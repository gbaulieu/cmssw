#include "RecoHI/HiTracking/interface/HIProtoTrackFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/Common/interface/DetSetAlgorithm.h"

#include "DataFormats/Common/interface/DetSetVector.h"    
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

using namespace std;
using namespace edm;

/*****************************************************************************/
HIProtoTrackFilter::HIProtoTrackFilter (const edm::ParameterSet& ps, const edm::EventSetup& es) :
theTIPMax( ps.getParameter<double>("tipMax") ),
theChi2Max( ps.getParameter<double>("chi2") ),
thePtMin( ps.getParameter<double>("ptMin") ),
doVariablePtMin( ps.getParameter<bool>("doVariablePtMin") ),
theBeamSpotTag( ps.getParameter<InputTag>("beamSpot")),
theSiPixelRecHits( ps.getParameter<InputTag>("siPixelRecHits")),
theBeamSpot(0),
theTrackerTopology(0),
theVariablePtMin(0)
{ 
  edm::ESHandle<TrackerTopology> httopo;
  es.get<IdealGeometryRecord>().get(httopo);
  theTrackerTopology = httopo.product();
}

/*****************************************************************************/
HIProtoTrackFilter::HIProtoTrackFilter (const edm::ParameterSet& ps) :
theTIPMax( ps.getParameter<double>("tipMax") ),
theChi2Max( ps.getParameter<double>("chi2") ),
thePtMin( ps.getParameter<double>("ptMin") ),
doVariablePtMin( ps.getParameter<bool>("doVariablePtMin") ),
theBeamSpotTag( ps.getParameter<InputTag>("beamSpot")),
theSiPixelRecHits( ps.getParameter<InputTag>("siPixelRecHits")),
theBeamSpot(0),
theVariablePtMin(0)
{
  throw cms::Exception("Configuration") << "HIProtoTrackFilter needs EventSetup, please use PixelTrackFilterWithESFactory to construct it";
}

/*****************************************************************************/
HIProtoTrackFilter::~HIProtoTrackFilter()
{ }

/*****************************************************************************/
bool HIProtoTrackFilter::operator() (const reco::Track* track,const PixelTrackFilter::Hits & recHits) const
{

  if (!track) return false; 

  float minpt = thePtMin;
  if(doVariablePtMin) minpt = theVariablePtMin;

  if (track->chi2() > theChi2Max || track->pt() < minpt) return false; 
  
  math::XYZPoint vtxPoint(0.0,0.0,0.0);
  
  if(theBeamSpot)
    vtxPoint = theBeamSpot->position();
  
  double d0=0.0;
  d0 = -1.*track->dxy(vtxPoint);
  
  if (theTIPMax>0 && fabs(d0)>theTIPMax) return false;
  
  return true;
}

/*****************************************************************************/
void HIProtoTrackFilter::update(edm::Event& ev, const edm::EventSetup& es)
{
  
  // Get the beam spot
  edm::Handle<reco::BeamSpot> bsHandle;
  ev.getByLabel( theBeamSpotTag, bsHandle);
  theBeamSpot = bsHandle.product();
  
  if(theBeamSpot) {
    LogInfo("HeavyIonVertexing") 
      << "[HIProtoTrackFilter] Proto track selection based on beamspot"
      << "\n   (x,y,z) = (" << theBeamSpot->x0() << "," << theBeamSpot->y0() << "," << theBeamSpot->z0() << ")";
  } else {
    LogError("HeavyIonVertexing") // this can be made a warning when operator() is fixed
      << "No beamspot found with tag '" << theBeamSpotTag << "'";
  }

  // Estimate multiplicity
  edm::Handle<SiPixelRecHitCollection> recHitColl;
  ev.getByLabel(theSiPixelRecHits, recHitColl);

  vector<const TrackingRecHit*> theChosenHits; 	 
  edmNew::copyDetSetRange(*recHitColl,theChosenHits, theTrackerTopology->pxbDetIdLayerComparator(1));
  float estMult = theChosenHits.size();
  
  theVariablePtMin=thePtMin;

  // parameterize ptMin such that a roughly constant number of selected prototracks passed are to vertexing
  float varPtCutoff = 1500; //cutoff for variable ptMin
  if(estMult < varPtCutoff) {
    theVariablePtMin = 0.075;
    if(estMult > 0) theVariablePtMin = (13. - (varPtCutoff/estMult) )/12.; 
    if(theVariablePtMin<0.075) theVariablePtMin = 0.075; // don't lower the cut past 75 MeV
  }
  
  LogTrace("heavyIonHLTVertexing")<<"   [HIProtoTrackFilter: theVariablePtMin: " << theVariablePtMin << "]";
  
  
  return;
  
}
