/*! \class   TrackFindingAMProducer
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 *  \update by S.Viret 
 *  \date   2014, Feb 17
 *
 */

#ifndef TRACK_BUILDER_AM_H
#define TRACK_BUILDER_AM_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "L1Trigger/TrackFindingAM/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackFindingAM/interface/PatternFinder.h"
#include "L1Trigger/TrackFindingAM/interface/SectorTree.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#ifndef __APPLE__
BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer)
#endif

class TrackFindingAMProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TrackFindingAMProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TrackFindingAMProducer();

  private:

  /// Data members
  std::string                  nBKName;
  int                          nThresh;
  int                          nMissingHits;
  int                          nDebug;
  SectorTree                   m_st;
  PatternFinder                *m_pf;

  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > m_stoken;

  //  edm::InputTag                TTStubsInputTag;
  //  edm::InputTag                TTClustersInputTag;
  std::string                  TTPatternOutputTag;

  edm::ESHandle<TrackerTopology> tTopoHandle;
  edm::ESHandle<TrackerGeometry> tGeomHandle;

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

  bool inPattern(int j);

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
TrackFindingAMProducer::TrackFindingAMProducer( const edm::ParameterSet& iConfig )
{
  m_stoken = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter< edm::InputTag >( "TTInputStubs" ));
  TTPatternOutputTag = iConfig.getParameter< std::string >( "TTPatternName" );
  nBKName            = iConfig.getParameter< std::string >("inputBankFile");
  nThresh            = iConfig.getParameter< int >("threshold");
  nMissingHits       = iConfig.getParameter< int >("nbMissingHits");
  nDebug             = iConfig.getParameter< int >("debugMode");

  std::cout << "Loading pattern bank file : " << std::endl;
  std::cout << nBKName << std::endl;

  std::ifstream ifs(nBKName.c_str());

  //boost::archive::text_iarchive ia(ifs);
  boost::iostreams::filtering_stream<boost::iostreams::input> f;
  f.push(boost::iostreams::gzip_decompressor());
  try { 
    f.push(ifs);
    boost::archive::text_iarchive ia(f);
    ia >> m_st;
  }
  catch (boost::iostreams::gzip_error& e) {
    if(e.error()==4){//file is not compressed->read it without decompression
      std::ifstream new_ifs(nBKName.c_str());
      boost::archive::text_iarchive ia(new_ifs);
      ia >> m_st;
    }
  }  

  m_pf = new PatternFinder( nThresh, &m_st, "", "" );

  if(nMissingHits>-1)
  {
    m_pf->useMissingHitThreshold(nMissingHits);
  }

  produces< std::vector< TTTrack<  Ref_Phase2TrackerDigi_  > > >( TTPatternOutputTag );
}

/// Destructor
TrackFindingAMProducer::~TrackFindingAMProducer() {}

/// Begin run
void TrackFindingAMProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);
}

/// End run
void TrackFindingAMProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFindingAMProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();


  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  auto ttTracksForOutput = std::make_unique<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>();

  /// Get the Stubs/Cluster already stored
  edm::Handle< edmNew::DetSetVector< TTStub<  Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  iEvent.getByToken( m_stoken, TTStubHandle );

  /// STEP 0
  /// Prepare output

  ttTracksForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;

  /// STEP 1
  /// Loop over input stubs

  std::vector< Hit* > m_hits;
  for(unsigned int i=0;i<m_hits.size();i++) delete m_hits[i];
  m_hits.clear();

  unsigned int j = 0;
  std::map< unsigned int, edm::Ref< edmNew::DetSetVector< TTStub<  Ref_Phase2TrackerDigi_  > >, TTStub<  Ref_Phase2TrackerDigi_  > > > stubMap;

  //  std::cout << "Dealing with " << TTStubHandle->size() << " modules containing stub(s)" << std::endl;

  /// Go on only if there are stubs
  if ( TTStubHandle->size() > 0 )
  {
    /// Loop over TTStubs
    for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) 
    {
      DetId detid = (*gd)->geographicalId();
      if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
      if(!tTopo->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
      DetId stackDetid = tTopo->stack(detid); // Stub module detid

      if (TTStubHandle->find( stackDetid ) == TTStubHandle->end() ) continue;

      //      std::cout << "Here" << std::endl;

      /// Get the DetSets of the Stubs
      edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*TTStubHandle)[ stackDetid ];
      const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
      const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

      for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) 
      {
	/// Increment the counter
	j++;

	/// Make the Ref to be put in the Track
	edm::Ref< edmNew::DetSetVector< TTStub<  Ref_Phase2TrackerDigi_  > >, TTStub<  Ref_Phase2TrackerDigi_  > > tempStubRef = makeRefTo( TTStubHandle, stubIter );

	if (tempStubRef->getTriggerDisplacement()==500 && nDebug==1)  continue; // Don't pass FE innef (under test)

	//	std::cout << "Dealing with stub " << j << std::endl;

	stubMap.insert( std::make_pair( j, tempStubRef ) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint coords = tempStubRef->getClusterRef(0)->findAverageLocalCoordinates();

	LocalPoint clustlp   = topol->localPosition(coords);
	GlobalPoint posStub  =  theGeomDet->surface().toGlobal(clustlp);

	/// Find pixel pitch and topology related information
	int segment = floor( coords.y() ); // Segment of the bottom clust

	// Here we rearrange the number in order to be compatible with the AM emulator
	if ( detid.subdetId()==StripSubdetector::TOB )
	{
	  layer  = static_cast<int>(tTopo->layer(detid))+4;
	  ladder = static_cast<int>(tTopo->tobRod(detid));
	  module = static_cast<int>(tTopo->module(detid));
	}
	else if ( detid.subdetId()==StripSubdetector::TID )
	{	
	  layer  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
	  ladder = static_cast<int>(tTopo->tidRing(detid));
	  module = static_cast<int>(tTopo->module(detid));
	}

	module = CMSPatternLayer::getModuleCode(layer,module);

	// the stub is on the third Z position on the other side of the tracker -> out of range
	if ( module < 0 )  continue;

	ladder = CMSPatternLayer::getLadderCode(layer, ladder);
	segment = CMSPatternLayer::getSegmentCode(layer, ladder, segment);

	float x    = posStub.x();
	float y    = posStub.y();
	float z    = posStub.z();

	Hit* h = new Hit(layer,ladder, module, segment, coords.x(), j, -1, 0, 0, 0, 0, x, y, z, 0, 0, 0, tempStubRef->getTriggerDisplacement()-tempStubRef->getTriggerOffset());
	if(m_st.getSector(*h)!=NULL)
	  m_hits.push_back(h);
	else
	  delete(h);
      } /// End of loop over input stubs
    } /// End of loop over DetSetVector
  }

  if(nDebug==1)
    std::cout << "End of stub loop, collected " << j << " stubs" << std::endl;

  /// STEP 2
  /// PAssing the superstrips into the AM chip
  if(nDebug==1)
    m_pf->setVerboseMode(true); // display the supertrips of the event
  std::vector< Sector* > patternsSectors = m_pf->find(m_hits); // AM PR is done here....

  if(nDebug==1)
    std::cout<<"Looking "<<patternsSectors.size()<< " sectors " << endl;

  /// STEP 3
  /// Collect the info and store the track seed stuff

  std::vector< Hit* > hits;
  std::vector< char > layers_touched;

  std::vector< edm::Ref< edmNew::DetSetVector< TTStub<  Ref_Phase2TrackerDigi_  > >, TTStub<  Ref_Phase2TrackerDigi_  > > > tempVec;

  for ( unsigned int i = 0; i < patternsSectors.size(); i++ )
  {
    std::vector< GradedPattern* > pl = patternsSectors[i]->getPatternTree()->getLDPatterns();

    if ( pl.size() == 0 ) continue; // No patterns

    int secID = patternsSectors[i]->getOfficialID();
    if(nDebug==1)
      std::cout<<"Found "<<pl.size()<<" patterns in sector " << secID<<std::endl;
    //std::cout<<"containing "<<n_active<<" layers " << secID<<std::endl;
  
    //delete the GradedPattern objects

    for ( unsigned j = 0; j < pl.size(); j++ )
    {
      hits.clear();
      hits = pl[j]->getHits();

      //      std::cout << pl[j]->getOrderInChip() << std::endl;
      
      int n_lay = pl[j]->getNbLayers() - pl[j]->getNbFakeSuperstrips(); 
      bool is_there;

      // Look how many layers have been hit.
      layers_touched.clear();

      for(unsigned k = 0; k < hits.size(); k++ )
      {
	is_there = false;
	if (layers_touched.size()==0) 
	{
	  layers_touched.push_back( hits[k]->getLayer());
	}
	else
	{	 
	  for(unsigned l = 0; l < layers_touched.size(); l++ )
	  {
	    if (is_there) break;
	    if (hits[k]->getLayer() == layers_touched.at(l)) is_there=true;
	  } 

	  if (!is_there) layers_touched.push_back( hits[k]->getLayer());
        }
      }

      /// Create the Seed in the form of a Track and store it in the output
      tempVec.clear();

      for(unsigned k = 0; k < hits.size(); k++ )
        tempVec.push_back( stubMap[ hits[k]->getID() ] );


      TTTrack<  Ref_Phase2TrackerDigi_ > tempTrack( tempVec );
      tempTrack.setSector( 100*pl[j]->getOrderInChip()+secID );
      tempTrack.setWedge( n_lay - layers_touched.size() );
      tempTrack.setPOCA( GlobalPoint(0.,0.,0.),5);		
      ttTracksForOutput->push_back( tempTrack );

      delete pl[j];
    }

    //delete the Sectors
    delete patternsSectors[i];
  }

  /// Put in the event content
  iEvent.put( std::move(ttTracksForOutput), TTPatternOutputTag);
}
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFindingAMProducer);

#endif

