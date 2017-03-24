/*! \class   TrackFitTCProducer
 *
 *  \author S Viret / G Baulieu / G Galbit
 *  \date   2015, Mar 10
 *
 */


#ifndef TRACK_FITTER_AM_H
#define TRACK_FITTER_AM_H

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
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "L1Trigger/TrackFindingAM/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackFindingAM/interface/PatternFinder.h"
#include "L1Trigger/TrackFindingAM/interface/SectorTree.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

//#ifndef __APPLE__
//BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer)
//#endif

class TrackFitTCProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TrackFitTCProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TrackFitTCProducer();

  private:
  
  /// Data members
  unsigned int                 nSectors;
  unsigned int                 nWedges;
  std::string                  nBKName;
  int                          nThresh;

  std::string                  TTTrackOutputTag;
  std::string                  TTTrackBinaryOutputTag;

  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > m_stoken;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > m_ptoken;

  std::string                  TTPatternOutputTag;

  edm::ESHandle<TrackerTopology> tTopoHandle;
  edm::ESHandle<TrackerGeometry> tGeomHandle;

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
TrackFitTCProducer::TrackFitTCProducer( const edm::ParameterSet& iConfig )
{
  m_stoken = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter< edm::InputTag >( "TTInputStubs" ));
  m_ptoken = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter< edm::InputTag >( "TTInputPatterns" ));

  TTTrackOutputTag         = iConfig.getParameter< std::string >( "TTTrackName" );
  TTTrackBinaryOutputTag   = iConfig.getParameter< std::string >( "TTTrackBinaryName" );

  produces< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >( TTTrackOutputTag );
  produces< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >( TTTrackBinaryOutputTag );
}

/// Destructor
TrackFitTCProducer::~TrackFitTCProducer() {}

/// Begin run
void TrackFitTCProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);
}

/// End run
void TrackFitTCProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void TrackFitTCProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  auto ttTracksForOutput    = std::make_unique<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>();
  auto ttTracksBinForOutput = std::make_unique<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>();

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTPatternHandle;

  iEvent.getByToken( m_stoken, TTStubHandle );
  iEvent.getByToken( m_ptoken, TTPatternHandle );

  /// STEP 0
  /// Prepare output
  ttTracksForOutput->clear();
  ttTracksBinForOutput->clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;

  int nbLayers = 0;

  /// Loop over Patterns
  unsigned int tkCnt = 0;
  unsigned int j     = 0;

  std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubMap;
  

  TCBuilder* TCB  = new TCBuilder(nbLayers); // Floating point
  TCBuilder* TCBb = new TCBuilder(nbLayers); // Bit-wise
  TCBb->setHardwareEmulation(true);

  /// STEP 1
  /// Loop over patterns

  //  std::cout << "Start the loop over " << TTPatternHandle->size() << " pattern(s) in order to recover the stubs" << std::endl;

  std::vector<Hit*> m_hits;
  for(unsigned int i=0;i<m_hits.size();i++) delete m_hits[i];
  m_hits.clear();

  std::vector<Track*> tracks;
  for(unsigned int i=0;i<tracks.size();i++) delete tracks[i];
  tracks.clear();

  std::vector<Track*> tracksb;
  for(unsigned int i=0;i<tracksb.size();i++) delete tracksb[i];
  tracksb.clear();

  edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator stubIter;

  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterTTTrack;
  
  /// Go on only if there are Patterns from PixelDigis
  if ( TTPatternHandle->size() > 0 )
  {
    for ( iterTTTrack = TTPatternHandle->begin();
	  iterTTTrack != TTPatternHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > tempTrackPtr( TTPatternHandle, tkCnt++ );

      j = 0;

      m_hits.clear();
      tracks.clear();
      stubMap.clear();

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();
      nbLayers = tempTrackPtr->getWedge();

      // Get the stubs in the road

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      // Loop over stubs contained in the pattern to recover the info

      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	++j;

	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = trackStubs.at(i);

	stubMap.insert( std::make_pair( j, tempStubRef ) );

	DetId stackDetid = tempStubRef->getDetId();
	DetId detid = stackDetid;

	// DIRTY!! but need this loop to get back the geographicalId from the detid 
	// Is there a method in TrackerTopology.h for that????
	for (auto gd=tGeomHandle->dets().begin(); gd != tGeomHandle->dets().end(); gd++) 
	{
	    DetId detidg = (*gd)->geographicalId();
	    if(detidg.subdetId()!=StripSubdetector::TOB && detidg.subdetId()!=StripSubdetector::TID ) continue; // only run on OT
	    if(tTopoHandle->stack(detidg)!=stackDetid) continue; 

	    detid=detidg;
	    break;
	}

	const GeomDetUnit* det0 = tGeomHandle->idToDetUnit( detid );
	const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );

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
	  layer  = static_cast<int>(tTopoHandle->layer(detid))+4;
	  ladder = static_cast<int>(tTopoHandle->tobRod(detid));
	  module = static_cast<int>(tTopoHandle->module(detid));
	}
	else if ( detid.subdetId()==StripSubdetector::TID )
	{	
	  layer  = 10+static_cast<int>(tTopoHandle->tidWheel(detid))+abs(2-static_cast<int>(tTopoHandle->side(detid)))*7;
	  ladder = static_cast<int>(tTopoHandle->tidRing(detid));
	  module = static_cast<int>(tTopoHandle->module(detid));
	}

	
	/// Find the z-segment (back to 0/1)
	segment = CMSPatternLayer::getSegmentCode(layer, ladder,floor( coords.y() ));


	//	cout << layer << " / " << ladder << " / " << module << " / " << std::endl;

	int strip  =  coords.x();
	int tp     = -1;
	float eta  = 0;
	float phi0 = 0;
	float spt  = 0;
	float x    = posStub.x();
	float y    = posStub.y();
	float z    = posStub.z();
	float x0   = 0.;
	float y0   = 0.;
	float z0   = 0.;
	float ip   = sqrt(x0*x0+y0*y0);
	
	Hit* h = new Hit(layer,ladder, module, segment, strip, 
			 j, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0, tempStubRef->getTriggerDisplacement()-tempStubRef->getTriggerOffset());
	m_hits.push_back(h);

      } /// End of loop over track stubs

      TCB->setSectorID(tempTrackPtr->getSector()%100);
      TCB->fit(m_hits);
      TCBb->setSectorID(tempTrackPtr->getSector()%100);
      TCBb->fit(m_hits);

      tracks = TCB->getTracks();
      TCB->clean();
      tracksb = TCBb->getTracks();
      TCBb->clean();


      if (tracks.size()>1 || tracksb.size()>1) 
      {
	// Not a normal behaviour, TC builder should produce only 1 TC max per road
	continue;
      }
      
      // Store the tracks (no duplicate cleaning yet)
      //      cout<<"Found "<<tracks.size()<<" track"<<endl;

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > tempVec;

      for(unsigned int tt=0;tt<tracks.size();tt++)
      {	
	tempVec.clear();

	vector<int> stubs = tracks[tt]->getStubs();
	for(unsigned int sti=0;sti<stubs.size();sti++) tempVec.push_back( stubMap[ stubs[sti] ]);

	double pz = tracks[tt]->getCurve()/(tan(2*atan(exp(-tracks[tt]->getEta0()))));
	
	TTTrack< Ref_Phase2TrackerDigi_ > tempTrack( tempVec );
	GlobalPoint POCA(0.,0.,tracks[tt]->getZ0());
	GlobalVector mom(tracks[tt]->getCurve()*cos(tracks[tt]->getPhi0()),
			 tracks[tt]->getCurve()*sin(tracks[tt]->getPhi0()),
			 pz);
	
	//std::cout << tracks[tt]->getCurve()<<" / "<<tracks[tt]->getZ0() << " / " <<tracks[tt]->getEta0()<<" / "<<tracks[tt]->getPhi0()<< std::endl;
	
	tempTrack.setSector( seedSector );
	tempTrack.setWedge( tracks[tt]->getCharge() );
	tempTrack.setMomentum( mom , 5);
	tempTrack.setPOCA( POCA , 5);
	//std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << " / " << tempTrack.getPOCA().z() << std::endl;
	ttTracksForOutput->push_back( tempTrack );
	
	delete tracks[tt];
      }

      //      std::cout << "PT6" << std::endl;

      for(unsigned int tt=0;tt<tracksb.size();tt++)
      {	
	tempVec.clear();

	vector<int> stubs = tracksb[tt]->getStubs();
	for(unsigned int sti=0;sti<stubs.size();sti++) tempVec.push_back( stubMap[ stubs[sti] ]);

	double pz = tracksb[tt]->getCurve()/(tan(2*atan(exp(-tracksb[tt]->getEta0()))));
	
	TTTrack< Ref_Phase2TrackerDigi_ > tempTrack( tempVec );
	GlobalPoint POCA(0.,0.,tracksb[tt]->getZ0());
	GlobalVector mom(tracksb[tt]->getCurve()*cos(tracksb[tt]->getPhi0()),
			 tracksb[tt]->getCurve()*sin(tracksb[tt]->getPhi0()),
			 pz);
	
	//std::cout << tracks[tt]->getCurve()<<" / "<<tracks[tt]->getZ0() << " / " <<tracks[tt]->getEta0()<<" / "<<tracks[tt]->getPhi0()<< std::endl;
	
	tempTrack.setSector( seedSector );
	tempTrack.setWedge( -1 );
	tempTrack.setMomentum( mom , 5);
	tempTrack.setPOCA( POCA , 5);
	//std::cout << tracks[tt]->getZ0() << " / " << POCA.z() << " / " << tempTrack.getPOCA().z() << std::endl;
	ttTracksBinForOutput->push_back( tempTrack );
	
	delete tracksb[tt];
      }

      //      std::cout << "PT7" << std::endl;

    } // End of loop over patterns
      
    //    std::cout << "PT8" << std::endl;

    delete(TCB);    
    delete(TCBb);  

  }

  /// Put in the event content
  iEvent.put( std::move(ttTracksForOutput), TTTrackOutputTag);
  iEvent.put( std::move(ttTracksBinForOutput), TTTrackBinaryOutputTag);
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(TrackFitTCProducer);

#endif

