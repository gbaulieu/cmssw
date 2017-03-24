/*! \class   AMOutputMerger
 *
 *
 *  \update by S.Viret 
 *  \date   2014, Feb 17
 *
 */

#ifndef AM_MERGER_H
#define AM_MERGER_H

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
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

class AMOutputMerger : public edm::EDProducer
{
  public:
    /// Constructor
    explicit AMOutputMerger( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~AMOutputMerger();

  private:

  /// Data members
  std::string                   TTStubOutputTag;
  std::vector< edm::InputTag >  TTPatternsInputTags;
  std::string                   TTPatternOutputTag;
  std::vector<int>              stored_IDs;


  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > m_stoken;
  edm::EDGetTokenT< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > m_ctoken;
  std::vector<edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > > m_ptokens;

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
AMOutputMerger::AMOutputMerger( const edm::ParameterSet& iConfig )
{
  m_ctoken = consumes< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter< edm::InputTag >( "TTInputClusters" ));
  m_stoken = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > > >(iConfig.getParameter< edm::InputTag >( "TTInputStubs" ));

  TTPatternsInputTags = iConfig.getParameter< std::vector< edm::InputTag > >( "TTInputPatterns" );

  for ( auto iTag =  TTPatternsInputTags.begin(); iTag!=  TTPatternsInputTags.end(); iTag++ )
  {
    m_ptokens.push_back(consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_  > > >(*iTag));
  }

  TTStubOutputTag     = iConfig.getParameter< std::string >( "TTFiltStubsName" );
  TTPatternOutputTag  = iConfig.getParameter< std::string >( "TTPatternsName" );

  produces< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >( TTPatternOutputTag );
  produces< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >( TTStubOutputTag );
}

/// Destructor
AMOutputMerger::~AMOutputMerger() {}

/// Begin run
void AMOutputMerger::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references

  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);
}

/// End run
void AMOutputMerger::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void AMOutputMerger::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output

  /// Get the Stubs/Cluster already stored
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  iEvent.getByToken( m_stoken, TTStubHandle );

  edm::Handle< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > TTClusterHandle;
  iEvent.getByToken( m_ctoken, TTClusterHandle );

  // The container for filtered patterns / stubs 

  auto ttTracksForOutput = std::make_unique<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>();
  auto ttStubsForOutput  = std::make_unique<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>>();

  std::vector< edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > > TTPatternHandle;

  TTPatternHandle.clear();
  TTPatternHandle.resize(static_cast<int>(TTPatternsInputTags.size()));

  int m=0;

  for ( auto iTag =  m_ptokens.begin(); iTag!=  m_ptokens.end(); iTag++ )
  {
    iEvent.getByToken( *iTag, TTPatternHandle.at(m) );
    ++m;
  }

  //
  // First step is pretty simple, we just create a map of all the stubs 
  // contained in the event
  //
 
  unsigned int stub_n = 0;
  std::map< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > , unsigned int > stubMap;

  /// Loop over L1TkStubs
  for (auto gd=tGeomHandle->dets().begin(); gd != tGeomHandle->dets().end(); gd++) 
  {
    DetId detid = (*gd)->geographicalId();
    if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; 
    if(!tTopoHandle->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
    
    DetId stackDetid = tTopoHandle->stack(detid); // Stub module detid

    if (TTStubHandle->find( stackDetid ) == TTStubHandle->end() ) continue;

    /// Get the DetSets of the Clusters
    edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*TTStubHandle)[ stackDetid ];

    for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) 
    {
      /// Make the Ref to be put in the Track
      edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > tempStubRef = edmNew::makeRefTo( TTStubHandle, stubIter );
      ++stub_n;
      stubMap.insert( std::make_pair( tempStubRef, stub_n ) );
    }
  }
  
  //
  // In the second step, we merge all the patterns into a single container
  // because they are stored in vectors
  //
 
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterTTTrack;

  for ( unsigned j = 0; j < TTPatternsInputTags.size(); ++j )
  {
    edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTPatterns = TTPatternHandle.at(j);

    if ( TTPatterns->size() > 0 )
    {      
      for ( iterTTTrack = TTPatterns->begin();
	    iterTTTrack != TTPatterns->end();
	    ++iterTTTrack )
      {

	TTTrack< Ref_Phase2TrackerDigi_ > tempTTPatt(iterTTTrack->getStubRefs());
	
	tempTTPatt.setSector(iterTTTrack->getSector());
	tempTTPatt.setWedge(iterTTTrack->getWedge());
	tempTTPatt.setMomentum(iterTTTrack->getMomentum(5),5);
	tempTTPatt.setPOCA(iterTTTrack->getPOCA(5),5);
	tempTTPatt.setRInv(iterTTTrack->getRInv(5),5);
	tempTTPatt.setChi2(iterTTTrack->getChi2(5),5);

	ttTracksForOutput->push_back(tempTTPatt);
      }
    }
  }

  // Get the OrphanHandle of the accepted patterns
  
  edm::OrphanHandle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTPatternAcceptedHandle = iEvent.put( std::move(ttTracksForOutput), TTPatternOutputTag );
 
  
  //
  // Third step, we flag the stubs contained in the patterns stored
  // 
  //

  bool found;

  stored_IDs.clear();

  if ( TTPatternAcceptedHandle->size() > 0 )
  {
    /// Loop over Patterns
    unsigned int tkCnt = 0;

    for ( iterTTTrack = TTPatternAcceptedHandle->begin();
	  iterTTTrack != TTPatternAcceptedHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > tempTrackPtr( TTPatternAcceptedHandle, tkCnt++ );

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      // Loop over stubs contained in the pattern to recover the info

      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	found=false;

	for(unsigned int l = 0; l < stored_IDs.size(); ++l )
	{
	  if (found) continue;
	  if (stored_IDs.at(l)==int(stubMap[ trackStubs.at(i) ])) found=true;
	}

	if (!found) stored_IDs.push_back(stubMap[ trackStubs.at(i) ]);

      }
    }
  }
  
  //
  // Last step, we recreate the filtered stub container from there
  // 
  //

  unsigned int j2 = 0;

  for (auto gd=tGeomHandle->dets().begin(); gd != tGeomHandle->dets().end(); gd++) 
  {
    DetId detid = (*gd)->geographicalId();
    if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue; 
    if(!tTopoHandle->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
    
    DetId stackDetid = tTopoHandle->stack(detid); // Stub module detid

    if (TTStubHandle->find( stackDetid ) == TTStubHandle->end() ) continue;

    /// Get the DetSets of the Clusters
    edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*TTStubHandle)[ stackDetid ];

    /// Create the vector of stubs to be passed to the FastFiller
    std::vector< TTStub< Ref_Phase2TrackerDigi_ > > *tempOutput = new std::vector< TTStub< Ref_Phase2TrackerDigi_ > >();
    tempOutput->clear();

    for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) 
    {
      ++j2;

      if (!AMOutputMerger::inPattern(j2)) continue;

      TTStub< Ref_Phase2TrackerDigi_ > tempTTStub( stubIter->getDetId() );

      tempTTStub.addClusterRef(stubIter->getClusterRef(0));
      tempTTStub.addClusterRef(stubIter->getClusterRef(1));
      tempTTStub.setTriggerDisplacement( stubIter->getTriggerDisplacement() );
      tempTTStub.setTriggerOffset( stubIter->getTriggerOffset() );
      tempOutput->push_back( tempTTStub );
    }

    /// Create the FastFiller

    if ( tempOutput->size() > 0 )
    {
      typename edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::FastFiller tempOutputFiller( *ttStubsForOutput, stackDetid );
      for ( unsigned int m = 0; m < tempOutput->size(); m++ )
      {
        tempOutputFiller.push_back( tempOutput->at(m) );
      }
      if ( tempOutputFiller.empty() )
        tempOutputFiller.abort();
    }
  }
  
  /// Put in the event content
  iEvent.put( std::move(ttStubsForOutput), TTStubOutputTag);  
}



bool AMOutputMerger::inPattern(int j)
{
  for(unsigned l = 0; l < stored_IDs.size(); ++l )
  {    
    if (stored_IDs.at(l)==j) return true;
  }

  return false;
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(AMOutputMerger);

#endif

