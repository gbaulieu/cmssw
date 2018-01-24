/*! \class   RoadCleaningProducer
 *
 *  \author S Viret / G Baulieu / G Galbit
 *  \date   2015, Mar 10
 *
 */


#ifndef ROAD_CLEANER_H
#define ROAD_CLEANER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

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

#include "TFile.h"
#include "TChain.h"

//#ifndef __APPLE__
//BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer)
//#endif

class RoadCleaningProducer : public edm::EDProducer
{
  public:
    /// Constructor
    explicit RoadCleaningProducer( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~RoadCleaningProducer();

  private:
  
  /// Data members
  unsigned int                 nSectors;
  unsigned int                 nWedges;
  int                          nThresh;
  std::string                  TTTrackOutputTag;

  TChain*         m_ptdata;
  std::string     m_ptFile;
  bool            m_filter;
  
  std::vector<std::vector <float> > pt_patt_full,pt_charge_full;
  std::vector<float> pt_patt,pt_charge;
  float m_pto;
  int   m_order;
  int   m_charge;
  int   m_towid;
  bool  m_tilted;


  int  limits[6][3];

  edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > m_stoken;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > m_ptoken;


  edm::ESHandle<TrackerTopology> tTopoHandle;
  edm::ESHandle<TrackerGeometry> tGeomHandle;

  bool stubFilter(int lay, int lad, double sw, double ppt, int charge);    

  /// Mandatory methods
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

}; /// Close class

/*! \brief   Implementation of methods
 */

/// Constructors
RoadCleaningProducer::RoadCleaningProducer( const edm::ParameterSet& iConfig )
{
  m_stoken = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter< edm::InputTag >( "TTInputStubs" ));
  m_ptoken = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter< edm::InputTag >( "TTInputPatterns" ));

  TTTrackOutputTag         = iConfig.getParameter< std::string >( "TTTrackName" );
  m_filter                 = iConfig.getParameter< bool >( "applyCuts" );
  m_ptFile                 = iConfig.getParameter< std::string >( "ptRoadFile" );
  m_tilted                 = !(iConfig.getParameter<bool>("flatBarrel"));

  produces< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >( TTTrackOutputTag );

  int n_tilted_rings[6];
  int n_flat_rings[6];

  for (int i=0; i < 6; ++i) n_tilted_rings[i]=0;
  for (int i=0; i < 6; ++i) n_flat_rings[i]=0;

  if (m_tilted)
  {
    n_tilted_rings[0]=11;
    n_tilted_rings[1]=12;
    n_tilted_rings[2]=12;
    n_flat_rings[0]=7;
    n_flat_rings[1]=11;
    n_flat_rings[2]=15;
  }

  for (int i=0; i < 6; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      limits[i][j]=0;

      if (n_tilted_rings[i]==0) continue;

      limits[i][j]=(j%2)*n_flat_rings[i]+(j>0)*n_tilted_rings[i];
    }
  }
}

/// Destructor
RoadCleaningProducer::~RoadCleaningProducer() {}

/// Begin run
void RoadCleaningProducer::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

  if (m_filter)
  {
    m_ptdata = new TChain("PatternData");
    m_ptdata->Add(m_ptFile.c_str());

    m_ptdata->SetBranchAddress("AveragePT",        &m_pto);
    m_ptdata->SetBranchAddress("order",            &m_order);
    m_ptdata->SetBranchAddress("charge",           &m_charge);
    m_ptdata->SetBranchAddress("sector",           &m_towid);

    pt_patt_full.clear();
    pt_charge_full.clear();
    
    pt_patt.clear();
    pt_charge.clear();

    for (int i=0;i<48;++i)
    {
      pt_patt_full.push_back(pt_patt);
      pt_charge_full.push_back(pt_charge);
    }

    for (unsigned int i=0;i<m_ptdata->GetEntries();++i)
    {
      m_ptdata->GetEntry(i);
      pt_patt_full.at(m_towid).push_back(0);
      pt_charge_full.at(m_towid).push_back(0);
    }
    
    for (unsigned int i=0;i<m_ptdata->GetEntries();++i)
    {
      m_ptdata->GetEntry(i);
      
      //      std::cout << pt_patt_full.at(m_towid).size() << " / " << m_towid << " / " << m_order << std::endl;
      
      pt_patt_full.at(m_towid).at(m_order)=m_pto;
      pt_charge_full.at(m_towid).at(m_order)=m_charge;    
    }
    
    std::cout << "PT info of banks has been stored... " << std::endl;
    std::cout << "This store contains " << m_ptdata->GetEntries() << " references..." << std::endl;
    std::cout << "Now filters the stubs from matched roads... " << std::endl;
  }
}

/// End run
void RoadCleaningProducer::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
void RoadCleaningProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  const TrackerTopology* const tTopo = tTopoHandle.product();
  const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  auto TTTracksForOutput = std::make_unique<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>();

  /// Get the Stubs already stored away
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTPatternHandle;

  iEvent.getByToken( m_stoken, TTStubHandle );
  iEvent.getByToken( m_ptoken, TTPatternHandle );

  /// STEP 0
  /// Prepare output
  TTTracksForOutput->clear();
  

  float barrel_sw_low[6] = {
    0.0, 0.0, 0.5, 1, 1.5, 2.0}; 


  float barrel_sw[6][5] = {
    { 1.7, 1.6, 1.6, 1.6, 1.6}, 
    { 1.6, 1.5, 1.5, 1.5, 1.4}, 
    { 2.2, 1.8, 1.7, 1.6, 1.5}, 
    { 3.3, 2.4, 2.2, 2.1, 2.2}, 
    { 4.1, 2.9, 3.7, 4, 3.2}, 
    { 4.8, 3.6, 4.8, 5.4, 4.4}
  }; 

  float disks_sw[15][5] = {
    { 1.3, 1.2, 1.1, 1.2, 1.1},{ 1.3, 1.3, 1.2, 1.2, 1.2},{ 1.5, 1.4, 1.3, 1.3, 1.3},{ 1.7, 1.5, 1.4, 1.4, 1.4},{ 1.8, 1.6, 1.5, 1.6, 1.4},{ 2, 1.7, 1.5, 1.5, 1.5},{ 2.1, 1.8, 1.6, 1.5, 1.5},{ 2.3, 1.9, 1.7, 1.6, 1.5},{ 2.5, 2, 1.7, 1.6, 1.5},{ 3.2, 2.5, 2, 1.8, 1.7},{ 2.9, 2.3, 2, 1.8, 1.7},{ 2.6, 2.1, 1.7, 1.5, 1.5},{ 2.7, 2, 1.6, 1.4, 1.4},{ 3.2, 2.3, 1.7, 1.8, 1.3},{ 3.6, 2.6, 1.8, 1.6, 1.3}
  };  

  int lay  = 0;
  int lad = 0;
  int nlaymiss = 0;
  int type;

  /// Loop over Patterns
  unsigned int tkCnt = 0;
  unsigned int j     = 0;

int inco;

  std::map< unsigned int , edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubMap;
  
  /// STEP 1
  /// Loop over patterns

  edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator stubIter;

  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterTTTrack;
  
  /// Go on only if there are Patterns from Phase2TrackerDigis
  if ( TTPatternHandle->size() > 0 )
  {
    for ( iterTTTrack = TTPatternHandle->begin();
	  iterTTTrack != TTPatternHandle->end();
	  ++iterTTTrack )
    {
      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > tempTrackPtr( TTPatternHandle, tkCnt++ );

      j = 0;
      inco=0;
      stubMap.clear();

      /// Get everything relevant
      unsigned int seedSector = tempTrackPtr->getSector();
      int realSector = static_cast<int>(seedSector)%100;
      nlaymiss = tempTrackPtr->getWedge();

      float ppt=0;
      float charge=0;

      if (m_filter)
      {
	ppt    = pt_patt_full.at(realSector).at(seedSector/100);
	charge = -pt_charge_full.at(realSector).at(seedSector/100);
      }      

      std::vector<int> layers_touched;
      std::vector<int> layers_touched_full;
      std::vector<int> filt_stubs;

      // Look how many layers have been hit.
      layers_touched.clear();
      layers_touched_full.clear();
      filt_stubs.clear();

      bool is_there;

      // Get the stubs in the road

      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > > > trackStubs = tempTrackPtr->getStubRefs();

      for(unsigned int i=0;i<trackStubs.size();i++)
      {
	++j;

	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = trackStubs.at(i);

	DetId stackDetid = tempStubRef->getDetId();
	DetId detid = stackDetid;

	stubMap.insert( std::make_pair( j, tempStubRef ) );

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
	const GeomDetUnit* det1 = tGeomHandle->idToDetUnit( detid+1 );

	const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
	const PixelTopology* top0 = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );

	/// Calculate average coordinates col/row for inner/outer Cluster
	/// These are already corrected for being at the center of each pixel
	MeasurementPoint coords = tempStubRef->getClusterRef(0)->findAverageLocalCoordinatesCentered();
	LocalPoint clustlp   = top0->localPosition(coords);
	GlobalPoint posStub  = pix0->surface().toGlobal(clustlp);
     
	double displStub     = tempStubRef->getTriggerDisplacement();
	double offsetStub    = tempStubRef->getTriggerOffset();

	int strip   = coords.x();
	int nstrips = top0->nrows();

	/// Get the Stack radius and z and displacements
	double R0 = det0->position().perp();
	double R1 = det1->position().perp();
	double Z0 = det0->position().z();
	double Z1 = det1->position().z();

	double DR = R1-R0;
	double DZ = Z1-Z0;

	double tilt  = atan2(DR,DZ);
	//double thick = sqrt(DR*DR+DZ*DZ)/(R0*sin(tilt)+Z0*cos(tilt));
	double thick = sqrt(DR*DR+DZ*DZ);

	// Here we rearrange the number in order to be compatible with the AM emulator
	if ( detid.subdetId()==StripSubdetector::TOB )
	{
	  type   = static_cast<int>(tTopo->tobSide(detid)); // Tilt-/Tilt+/Flat <-> 1/2/3
	  lay  = static_cast<int>(tTopo->layer(detid))+4;
	  lad = static_cast<int>(tTopo->tobRod(detid))-1;

	  if (type<3)
	  {
	    lad = static_cast<int>(tTopo->module(detid))-1;
	  }
	}
	else if ( detid.subdetId()==StripSubdetector::TID )
	{	
	  lay  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
	  lad = static_cast<int>(tTopo->tidRing(detid))-1;

	  // Special treatment for Tilted4021 geometry
	  if (static_cast<int>(tTopo->tidWheel(detid))>=3 && m_tilted) lad += 2;

	  type   = 0;
	}



	float x    = posStub.x();
	float y    = posStub.y();
	float z    = posStub.z();
	float r    = sqrt(x*x+y*y);

	// Get the real parallax correction
	double par_real=(strip-nstrips/2)*thick/(r*sin(tilt)+z*cos(tilt));

	double par_cor = par_real-offsetStub; // The correction to the correction

	//if (lay>=18) par_cor = -par_real-offsetStub;

	// Now compute the corrected stub width 

	double sw = displStub-offsetStub-par_cor;

	if (lay>=11 && lay<18)
	{
	  sw = -displStub+offsetStub+par_cor;
	}

	if (fabs(par_cor)>0.55)
	{
	  std::cout << lay << " / " << lad << " / " << strip << " / " 
		    << tilt << " / " << thick << " / " << par_real << " / " << offsetStub << " / " 
		    << par_cor << std::endl;
	}
      
	is_there = false;
	if (layers_touched_full.size()==0) 
	{
	  layers_touched_full.push_back( lay);
	}
	else
	{	 
	  for(unsigned l = 0; l < layers_touched_full.size(); l++ )
	  {
	    if (is_there) break;
	    if (lay == layers_touched_full.at(l)) is_there=true;
	  } 

	  if (!is_there) layers_touched_full.push_back( lay );
        }
	
	if (m_filter)
	{	  	  
	  int ptbin = std::min(int(ppt)/5,4);
    
	  if (lay>=11)
	  {
	    if (std::fabs(sw)>disks_sw[lad][ptbin]) continue;
	    //if (std::fabs(sw)*ppt>disks_pr[lad][ptbin]) continue;	    
	    //if (charge*sw<=0 && std::fabs(sw)>disks_ch[lad][ptbin]) continue;

	    if (charge*sw<=0 && lad<9 && std::fabs(sw)>1.3) continue;
	    if (charge*sw<=0 && lad>=9 && std::fabs(sw)>1.5) continue;
	  }
	  else
	  {
	    if (ptbin==0 && std::fabs(sw)<barrel_sw_low[lay-5]) continue;	 
	    if (std::fabs(sw)>barrel_sw[lay-5][ptbin]) continue;	    
	    //if (std::fabs(sw)*ppt>barrel_pr[lay-5][ptbin]) continue;	  
	    //if (charge*sw<=0 && std::fabs(sw)>barrel_ch[lay-5][ptbin]) continue;

	    if (charge*sw<=0 && lay>7 && std::fabs(sw)>2.5) continue;
	    if (charge*sw<=0 && lad<=7 && std::fabs(sw)>1.6) continue;
	  }

	  //if (!RoadCleaningProducer::stubFilter(lay-5, lad, sw, ppt, charge)) continue;

	  if (ppt<=10 && charge*sw<=0) ++inco; 
	}
      
	filt_stubs.push_back(j);

	is_there = false;
	if (layers_touched.size()==0) 
	{
	  layers_touched.push_back( lay);
	}
	else
	{	 
	  for(unsigned l = 0; l < layers_touched.size(); l++ )
	  {
	    if (is_there) break;
	    if (lay == layers_touched.at(l)) is_there=true;
	  } 

	  if (!is_there) layers_touched.push_back( lay );
        }
      } /// End of loop over track stubs
     

      int missed = layers_touched_full.size()-layers_touched.size()+nlaymiss;

      // We never accept more than 1 layer missing wrt the original road size.

      if (inco>2) continue;
      if (missed>1) continue;
      
      std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > tempVec;

      tempVec.clear();

      for(unsigned int sti=0;sti<filt_stubs.size();sti++) tempVec.push_back( stubMap[ filt_stubs[sti] ]);
	
      if (tempVec.size()>15) continue;

      TTTrack< Ref_Phase2TrackerDigi_ > tempTrack( tempVec );
	
      tempTrack.setSector( seedSector );
      tempTrack.setWedge( missed );
      tempTrack.setPOCA( GlobalPoint(0.,0.,0.),5);
      TTTracksForOutput->push_back( tempTrack );
    } // End of loop over patterns
  }

  /// Put in the event content
  iEvent.put( std::move(TTTracksForOutput), TTTrackOutputTag);
}


/// This is where the stub filtering is done
bool RoadCleaningProducer::stubFilter(int lay, int lad, double sw, double ppt, int charge)
{ 
  /*
  // Basic stub width cuts for high pt tracks
               
  // Barrel
  if (lay==0 && fabs(sw)>1.6 && ppt>20) return false;
  if (lay==1 && fabs(sw)>1.2 && ppt>20) return false;
  if (lay==2 && fabs(sw)>1.5 && ppt>20) return false;
  if (lay==3 && fabs(sw)>3.0 && ppt>10) return false;
  if (lay==4 && fabs(sw)>1.8 && ppt>10) return false;
  if (lay==5 && fabs(sw)>2.2 && ppt>10) return false;
           
  // Disks
  if (lay>5 && lad==0 && fabs(sw)>0.8 && ppt>20) return false;
  if (lay>5 && lad==1 && fabs(sw)>1.0 && ppt>20) return false;
  if (lay>5 && lad==2 && fabs(sw)>1.1 && ppt>20) return false;
  if (lay>5 && lad==3 && fabs(sw)>1.2 && ppt>20) return false;
  if (lay>5 && lad==4 && fabs(sw)>1.5 && ppt>20) return false;
  if (lay>5 && lad==5 && fabs(sw)>1.2 && ppt>20) return false;
  if (lay>5 && lad==6 && fabs(sw)>1.2 && ppt>20) return false;
  if (lay>5 && lad==7 && fabs(sw)>1.2 && ppt>20) return false;
  if (lay>5 && lad==8 && fabs(sw)>1.2 && ppt>20) return false;
  if (lay>5 && lad==9 && fabs(sw)>2.5 && ppt>10) return false;
  if (lay>5 && lad==10 && fabs(sw)>2.0 && ppt>10) return false;
  if (lay>5 && lad==11 && fabs(sw)>2.8 && ppt>10) return false;
  if (lay>5 && lad==12 && fabs(sw)>2.0 && ppt>10) return false;
  if (lay>5 && lad==13 && fabs(sw)>2.2 && ppt>10) return false;
  if (lay>5 && lad==14 && fabs(sw)>1.5 && ppt>10) return false;
  

  // Cuts checking the consistency of the stub bend with the road pt
  // Barrel
  if (lay==0 && fabs(sw)*ppt>14 && ppt<20) return false;
  if (lay==1 && fabs(sw)*ppt>13 && ppt<20) return false;
  if (lay==2 && fabs(sw)*ppt>14 && ppt<20) return false;
  if (lay==3 && fabs(sw)*ppt>17 && ppt<10) return false;
  if (lay==4 && fabs(sw)*ppt>20 && ppt<10) return false;
  if (lay==5 && fabs(sw)*ppt>20 && ppt<10) return false;
  
  // Disks
  if (lay>5 && lad==0 && fabs(sw)*ppt>6 && ppt<20) return false;
  if (lay>5 && lad==1 && fabs(sw)*ppt>8 && ppt<20) return false;
  if (lay>5 && lad==2 && fabs(sw)*ppt>8 && ppt<20) return false;
  if (lay>5 && lad==3 && fabs(sw)*ppt>7 && ppt<20) return false;
  if (lay>5 && lad==4 && fabs(sw)*ppt>8 && ppt<20) return false;
  if (lay>5 && lad==5 && fabs(sw)*ppt>8 && ppt<20) return false;
  if (lay>5 && lad==6 && fabs(sw)*ppt>11 && ppt<20) return false;
  if (lay>5 && lad==7 && fabs(sw)*ppt>9 && ppt<20) return false;
  if (lay>5 && lad==8 && fabs(sw)*ppt>10 && ppt<20) return false;
  if (lay>5 && lad==9 && fabs(sw)*ppt>15 && ppt<10) return false;
  if (lay>5 && lad==10 && fabs(sw)*ppt>12 && ppt<10) return false;
  if (lay>5 && lad==11 && fabs(sw)*ppt>12 && ppt<10) return false;
  if (lay>5 && lad==12 && fabs(sw)*ppt>12 && ppt<10) return false;
  if (lay>5 && lad==13 && fabs(sw)*ppt>16 && ppt<10) return false;
  if (lay>5 && lad==14 && fabs(sw)*ppt>15 && ppt<10) return false;
          

  // Consistency of the bend with the road charge
  // Barrel
  if (charge*sw<0 && fabs(sw)>=1. && (lay==0)) return false;
  if (charge*sw<0 && fabs(sw)>=0.8 && (lay==1)) return false;
  if (charge*sw<0 && fabs(sw)>=0.7 && (lay==2)) return false;
  if (charge*sw<0 && fabs(sw)>=0.8 && (lay>2 && lay<=5)) return false;

  // Disk
  if (charge*sw<0 && fabs(sw)>=0.8 && (lay>5 && lad==9)) return false;
  if (charge*sw<0 && fabs(sw)>=0.6 && (lay>5 && lad==10)) return false;
  if (charge*sw<0 && fabs(sw)>=0.7 && (lay>5 && lad==11)) return false;
  if (charge*sw<0 && fabs(sw)>=0.6 && (lay>5 && lad==12)) return false;
  if (charge*sw<0 && fabs(sw)>=0.8 && (lay>5 && lad==13)) return false;
  if (charge*sw<0 && fabs(sw)>=0.6 && (lay>5 && lad==14)) return false;

  if (charge*sw<0 && fabs(sw)>=0.5 && (lay>5 && lad==8)) return false;
  if (charge*sw<0 && fabs(sw)>=0.6 && (lay>5 && lad==7)) return false;
  if (charge*sw<0 && fabs(sw)>=0.7 && (lay>5 && lad==6)) return false;
  if (charge*sw<0 && fabs(sw)>=0.7 && (lay>5 && lad==5)) return false;
  if (charge*sw<0 && fabs(sw)>=0.8 && (lay>5 && lad==4)) return false;
  if (charge*sw<0 && fabs(sw)>=0.7 && (lay>5 && lad==3)) return false;
  if (charge*sw<0 && fabs(sw)>=0.9 && (lay>5 && lad==2)) return false;
  if (charge*sw<0 && fabs(sw)>=0.9 && (lay>5 && lad==1)) return false;
  if (charge*sw<0 && fabs(sw)>=0.9 && (lay>5 && lad==8)) return false;
            
  //if (lay==0 && fabs(sw)>1.6) return false;
  if (lay==1 && fabs(sw)>1.7) return false;
  if (lay==2 && fabs(sw)>2.4) return false;
  if (lay==3 && fabs(sw)>3.5) return false;
  if (lay==4 && fabs(sw)>4.4) return false;
  if (lay==5 && fabs(sw)>4.8) return false;
                
  if (lay>5 && lad==0 && fabs(sw)>1.1) return false;
  if (lay>5 && lad==1 && fabs(sw)>1.2) return false;
  
  if (lay>5)
  {
    if (lad==2)
    {
      if ((lay==6 || lay==13) && fabs(sw)>1.4) return false;
      if ((lay==7 || lay==14) && fabs(sw)>1.2) return false;
      if ((lay==8 || lay==15) && fabs(sw)>1.) return false;
    }

    if (lad==3)
    {
      if ((lay==6 || lay==13) && fabs(sw)>1.5) return false;
      if ((lay==7 || lay==14) && fabs(sw)>1.4) return false;
      if ((lay==8 || lay==15) && fabs(sw)>1.3) return false;
    }
                
    if (lad==4)
    {
      if ((lay==6 || lay==13) && fabs(sw)>1.8) return false;
      if ((lay==7 || lay==14) && fabs(sw)>1.6) return false;
      if ((lay==8 || lay==15) && fabs(sw)>1.3) return false;
      if ((lay==9 || lay==16) && fabs(sw)>1.3) return false;
    }
                
    if (lad==5)
    {
      if ((lay==6 || lay==13) && fabs(sw)>2.) return false;
      if ((lay==7 || lay==14) && fabs(sw)>1.6) return false;
      if ((lay==8 || lay==15) && fabs(sw)>1.5) return false;
      if ((lay==9 || lay==16) && fabs(sw)>1.5) return false;
    }
                
    if (lad==6)
    {
      if ((lay==6  || lay==13) && fabs(sw)>2.) return false;
      if ((lay==7  || lay==14) && fabs(sw)>1.8) return false;
      if ((lay==8  || lay==15) && fabs(sw)>1.5) return false;
      if ((lay==9  || lay==16) && fabs(sw)>1.4) return false;
      if ((lay==10 || lay==17) && fabs(sw)>1.3) return false;
    }
    
    if (lad==7)
    {
      if ((lay==6  || lay==13) && fabs(sw)>2.1) return false;
      if ((lay==7  || lay==14) && fabs(sw)>1.8) return false;
      if ((lay==8  || lay==15) && fabs(sw)>1.8) return false;
      if ((lay==9  || lay==16) && fabs(sw)>1.7) return false;
      if ((lay==10 || lay==17) && fabs(sw)>1.6) return false;
    }
                
    if (lad==8)
    {
      if ((lay==6  || lay==13) && fabs(sw)>2.2) return false;
      if ((lay==7  || lay==14) && fabs(sw)>2.2) return false;
      if ((lay==8  || lay==15) && fabs(sw)>1.8) return false;
      if ((lay==9  || lay==16) && fabs(sw)>1.6) return false;
      if ((lay==10 || lay==17) && fabs(sw)>1.6) return false;
    }
                
    if (lad==9)
    {
      if ((lay==6  || lay==13) && fabs(sw)>3.3) return false;
      if ((lay==7  || lay==14) && fabs(sw)>2.8) return false;
      if ((lay==8  || lay==15) && fabs(sw)>2.4) return false;
      if ((lay==9  || lay==16) && fabs(sw)>2.1) return false;
      if ((lay==10 || lay==17) && fabs(sw)>1.7) return false;
    }
                
    if (lad==10)
    {
      if ((lay==6  || lay==13) && fabs(sw)>2.3) return false;
      if ((lay==7  || lay==14) && fabs(sw)>1.9) return false;
      if ((lay==8  || lay==15) && fabs(sw)>2.8) return false;
      if ((lay==9  || lay==16) && fabs(sw)>2.6) return false;
      if ((lay==10 || lay==17) && fabs(sw)>2.4) return false;
    }
    
    if (lad==11)
    {
      if ((lay==6  || lay==13) && fabs(sw)>2.5) return false;
      if ((lay==7  || lay==14) && fabs(sw)>2.2) return false;
      if ((lay==8  || lay==15) && fabs(sw)>2.0) return false;
      if ((lay==9  || lay==16) && fabs(sw)>1.5) return false;
      if ((lay==10 || lay==17) && fabs(sw)>2.6) return false;
    }
                
    if (lad==12)
    {
      if ((lay==6  || lay==13) && fabs(sw)>2.8) return false;
      if ((lay==7  || lay==14) && fabs(sw)>2.5) return false;
      if ((lay==8  || lay==15) && fabs(sw)>2.2) return false;
      if ((lay==9  || lay==16) && fabs(sw)>1.9) return false;
      if ((lay==10 || lay==17) && fabs(sw)>1.5) return false;
    }
                
    if (lad==13)
    {
      if ((lay==6  || lay==13) && fabs(sw)>3.2) return false;
      if ((lay==7  || lay==14) && fabs(sw)>2.7) return false;
      if ((lay==8  || lay==15) && fabs(sw)>2.5) return false;
      if ((lay==9  || lay==16) && fabs(sw)>2.2) return false;
      if ((lay==10 || lay==17) && fabs(sw)>1.6) return false;
    }
                
    if (lad==14)
    {
      if ((lay==6  || lay==13) && fabs(sw)>3.7) return false;
      if ((lay==7  || lay==14) && fabs(sw)>3.7) return false;
      if ((lay==8  || lay==15) && fabs(sw)>3.3) return false;
      if ((lay==9  || lay==16) && fabs(sw)>2.8) return false;
      if ((lay==10 || lay==17) && fabs(sw)>2.2) return false;
    }
  }
  */

  // Reoptimized
       
  if (lay==0 && fabs(sw)>1.7) return false;
  if (lay==1 && fabs(sw)>1.7) return false;
  if (lay==2 && fabs(sw)>2.3) return false;
  if (lay==3 && fabs(sw)>3.3) return false;
  if (lay==4 && fabs(sw)>4.1) return false;
  if (lay==5 && fabs(sw)>4.9) return false;

  if (charge==0)
  {
    if (lay==0 && fabs(sw)>1.5) return false;
    if (lay==1 && fabs(sw)>1.3) return false;
    if (lay==2 && fabs(sw)>1.2) return false;
    if (lay==3 && fabs(sw)>1.2) return false;
    if (lay==4 && fabs(sw)>1.2) return false;
    if (lay==5 && fabs(sw)>1.5) return false;
  }

        
  // Barrel
  if (lay==0 && fabs(sw)>1.5 && ppt>20) return false;
  if (lay==1 && fabs(sw)>1.4 && ppt>20) return false;
  if (lay==2 && fabs(sw)>1.4 && ppt>20) return false;
  if (lay==3 && fabs(sw)>2.0 && ppt>10) return false;
  if (lay==4 && fabs(sw)>3.0 && ppt>10) return false;
  if (lay==5 && fabs(sw)>4.1 && ppt>10) return false;
                 



  // Cuts checking the consistency of the stub bend with the road pt
  // Barrel
  if (lay==0 && fabs(sw)*ppt>12 && ppt<20) return false;
  if (lay==1 && fabs(sw)*ppt>10 && ppt<20) return false;
  if (lay==2 && fabs(sw)*ppt>11 && ppt<20) return false;
  if (lay==3 && fabs(sw)*ppt>15 && ppt<20) return false;
  if (lay==4 && fabs(sw)*ppt>18 && ppt<20) return false;
  if (lay==5 && fabs(sw)*ppt>21 && ppt<20) return false;
  

  // Consistency of the bend with the road charge
  if (charge*sw<0 && fabs(sw)>1.7) return false;
  if (lay>5 && charge==0 && fabs(sw)>1.7) return false; 

  if ((lay==6 || lay==13))
  {
    if (lad==0)
    {
      if (fabs(sw)>1.3) return false;
      if (fabs(sw)>1.1 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==1)
    {
      if (fabs(sw)>1.3) return false;
      if (fabs(sw)>1.1 && ppt>20)    return false;
      if (fabs(sw)*ppt>11 && ppt<20) return false;
    }

    if (lad==2)
    {
      if (fabs(sw)>1.5) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>11 && ppt<20) return false;
    }

    if (lad==3)
    {
      if (fabs(sw)>1.8) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==4)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.8 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==5)
    {
      if (fabs(sw)>2.1) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==6)
    {
      if (fabs(sw)>2.2) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==7)
    {
      if (fabs(sw)>2.5) return false;
      if (fabs(sw)>1.3 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==8)
    {
      if (fabs(sw)>2.8) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==9)
    {
      if (fabs(sw)>3.5) return false;
      if (fabs(sw)>1.9 && ppt>20)    return false;
      if (fabs(sw)*ppt>18 && ppt<20) return false;
    }

    if (lad==10)
    {
      if (fabs(sw)>2.3) return false;
      if (fabs(sw)>2.1 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==11)
    {
      if (fabs(sw)>2.6) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==12)
    {
      if (fabs(sw)>2.9) return false;
      if (fabs(sw)>1.9 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==13)
    {
      if (fabs(sw)>3.5) return false;
      if (fabs(sw)>1.8 && ppt>20)    return false;
      if (fabs(sw)*ppt>17 && ppt<20) return false;
    }

    if (lad==14)
    {
      if (fabs(sw)>3.8) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>17 && ppt<20) return false;
    }
  }

  if ((lay==7 || lay==14))
  {
    if (lad<1) return false;

    if (lad==1)
    {
      if (fabs(sw)>1.3) return false;
      if (fabs(sw)>1.1 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==2)
    {
      if (fabs(sw)>1.4) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==3)
    {
      if (fabs(sw)>1.6) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>11 && ppt<20) return false;
    }

    if (lad==4)
    {
      if (fabs(sw)>1.8) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==5)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==6)
    {
      if (fabs(sw)>2.) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==7)
    {
      if (fabs(sw)>2.1) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==8)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==9)
    {
      if (fabs(sw)>3.) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>16 && ppt<20) return false;
    }

    if (lad==10)
    {
      if (fabs(sw)>2.) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==11)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==12)
    {
      if (fabs(sw)>2.6) return false;
      if (fabs(sw)>2 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==13)
    {
      if (fabs(sw)>3.) return false;
      if (fabs(sw)>1.2 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==14)
    {
      if (fabs(sw)>3.4) return false;
      if (fabs(sw)>1.8 && ppt>20)    return false;
      if (fabs(sw)*ppt>16 && ppt<20) return false;
    }
  }

  if ((lay==8 || lay==15))
  {
    if (lad<2) return false;

    if (lad==2)
    {
      if (fabs(sw)>1.3) return false;
      if (fabs(sw)>1.3 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==3)
    {
      if (fabs(sw)>1.6) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==4)
    {
      if (fabs(sw)>1.6) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==5)
    {
      if (fabs(sw)>1.8) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==6)
    {
      if (fabs(sw)>1.8) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==7)
    {
      if (fabs(sw)>2) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==8)
    {
      if (fabs(sw)>2.2) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==9)
    {
      if (fabs(sw)>2.8) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>16 && ppt<20) return false;
    }

    if (lad==10)
    {
      if (fabs(sw)>3.1) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>17 && ppt<20) return false;
    }

    if (lad==11)
    {
      if (fabs(sw)>2.1) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==12)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==13)
    {
      if (fabs(sw)>2.6) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==14)
    {
      if (fabs(sw)>2.9) return false;
      if (fabs(sw)>1.3 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }
  }

  if ((lay==9 || lay==16))
  { 
    if (lad<3) return false;

    if (lad==3)
    {
      if (fabs(sw)>0.6) return false;
      if (fabs(sw)>0.3 && ppt>20)    return false;
      if (fabs(sw)*ppt>4 && ppt<20) return false;
    }

    if (lad==4)
    {
      if (fabs(sw)>1.5) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==5)
    {
      if (fabs(sw)>1.6) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==6)
    {
      if (fabs(sw)>1.6) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==7)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==8)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==9)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==10)
    {
      if (fabs(sw)>2.8) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>17 && ppt<20) return false;
    }

    if (lad==11)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==12)
    {
      if (fabs(sw)>2.1) return false;
      if (fabs(sw)>1.4 && ppt>20)    return false;
      if (fabs(sw)*ppt>12 && ppt<20) return false;
    }

    if (lad==13)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.8 && ppt>20)    return false;
      if (fabs(sw)*ppt>14 && ppt<20) return false;
    }

    if (lad==14)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.1 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }
  }

  if ((lay==10 || lay==17))
  {
    if (lad<5) return false;

    if (lad==5)
    {
      if (fabs(sw)>0.5) return false;
      if (fabs(sw)>0.5 && ppt>20)    return false;
      if (fabs(sw)*ppt>5 && ppt<20) return false;
    }

    if (lad==6)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.6 && ppt>20)    return false;
      if (fabs(sw)*ppt>16 && ppt<20) return false;
    }

    if (lad==7)
    {
      if (fabs(sw)>1.7) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==8)
    {
      if (fabs(sw)>1.7) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==9)
    {
      if (fabs(sw)>2.2) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==10)
    {
      if (fabs(sw)>2.4) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>15 && ppt<20) return false;
    }

    if (lad==11)
    {
      if (fabs(sw)>2.9) return false;
      if (fabs(sw)>1.7 && ppt>20)    return false;
      if (fabs(sw)*ppt>17 && ppt<20) return false;
    }

    if (lad==12)
    {
      if (fabs(sw)>1.9) return false;
      if (fabs(sw)>1.3 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==13)
    {
      if (fabs(sw)>2.3) return false;
      if (fabs(sw)>1.2 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }

    if (lad==14)
    {
      if (fabs(sw)>2.3) return false;
      if (fabs(sw)>1.2 && ppt>20)    return false;
      if (fabs(sw)*ppt>13 && ppt<20) return false;
    }
  }


  return true;
}

// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(RoadCleaningProducer);

#endif

