//
//
// Package:    Phase2TrackerDigitizer
// Class:      Phase2TrackerDigitizer
// 
// *\class SiPhase2TrackerDigitizer Phase2TrackerDigitizer.cc SimTracker/SiPhase2Digitizer/src/Phase2TrackerDigitizer.cc
//
// Description: <one line class summary>
//
// Implementation:
//     <Notes on implementation>
//

#include <memory>
#include <set>
#include <iostream>
#include "SimTracker/SiPhase2Digitizer/interface/Phase2TrackerDigitizer.h"
#include "SimTracker/SiPhase2Digitizer/interface/Phase2TrackerDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/interface/SSDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/interface/PSSDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/interface/PSPDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/interface/PixelDigitizerAlgorithm.h"
#include "SimTracker/SiPhase2Digitizer/interface/DigitizerUtility.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"

// Random Number
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CLHEP/Random/RandomEngine.h"

namespace cms
{
  const std::string Phase2TrackerDigitizer::InnerPixel = "P2Pixel"; 
  const std::string Phase2TrackerDigitizer::PixelinPS  = "PSP";
  const std::string Phase2TrackerDigitizer::StripinPS  = "PSS";
  const std::string Phase2TrackerDigitizer::TwoStrip   = "SS";   

  Phase2TrackerDigitizer::Phase2TrackerDigitizer(const edm::ParameterSet& iConfig, edm::EDProducer& mixMod):
    first_(true),
    hitsProducer_(iConfig.getParameter<std::string>("hitsProducer")),
    trackerContainers_(iConfig.getParameter<std::vector<std::string> >("ROUList")),
    geometryType_(iConfig.getParameter<std::string>("GeometryType"))
  {
    edm::LogInfo("Phase2TrackerDigitizer") << "Initialize Digitizer Algorithms";
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
        << "Phase2TrackerDigitizer requires the RandomNumberGeneratorService\n"
           "which is not present in the configuration file.  You must add the service\n"
           "in the configuration file or remove the modules that require it.";
    }
    rndEngine_ = &(rng->getEngine());

    // one type of Digi and DigiSimLink suffices 
    // changes in future: InnerPixel -> Tracker
    const std::string alias1("simSiPixelDigis"); 
    mixMod.produces<edm::DetSetVector<PixelDigi> >("Pixel").setBranchAlias(alias1);
    mixMod.produces<edm::DetSetVector<PixelDigiSimLink> >("Pixel").setBranchAlias(alias1);

    const std::string alias2("simSiTrackerDigis"); 
    mixMod.produces<edm::DetSetVector<Phase2TrackerDigi> >("Tracker").setBranchAlias(alias2);
    mixMod.produces<edm::DetSetVector<PixelDigiSimLink> >("Tracker").setBranchAlias(alias2);

    // creating algorithm objects and pushing them into the map
    algomap_[InnerPixel] = std::unique_ptr<Phase2TrackerDigitizerAlgorithm>(new PixelDigitizerAlgorithm(iConfig, (*rndEngine_)));
    algomap_[PixelinPS]  = std::unique_ptr<Phase2TrackerDigitizerAlgorithm>(new PSPDigitizerAlgorithm(iConfig, (*rndEngine_)));
    algomap_[StripinPS]  = std::unique_ptr<Phase2TrackerDigitizerAlgorithm>(new PSSDigitizerAlgorithm(iConfig, (*rndEngine_)));
    algomap_[TwoStrip]   = std::unique_ptr<Phase2TrackerDigitizerAlgorithm>(new SSDigitizerAlgorithm(iConfig, (*rndEngine_)));
  }
  Phase2TrackerDigitizer::~Phase2TrackerDigitizer() {  
    edm::LogInfo("Phase2TrackerDigitizer") << "Destroying the Digitizer";
  }
  void
  Phase2TrackerDigitizer::accumulatePixelHits(edm::Handle<std::vector<PSimHit> > hSimHits,
				       size_t globalSimHitIndex,const unsigned int tofBin) {
    if (hSimHits.isValid()) {
      std::set<unsigned int> detIds;
      std::vector<PSimHit> const& simHits = *hSimHits.product();
      for (auto it = simHits.begin(), itEnd = simHits.end(); it != itEnd; ++it, ++globalSimHitIndex) {
	unsigned int detId_raw = (*it).detUnitId();
        if (detectorUnits_.find(detId_raw) == detectorUnits_.end()) continue;
	if (detIds.insert(detId_raw).second) {
	  // The insert succeeded, so this detector element has not yet been processed.
	  const std::string algotype = getAlgoType(detId_raw);
	  Phase2TrackerGeomDetUnit* phase2det = detectorUnits_[detId_raw];
	  // access to magnetic field in global coordinates
	  GlobalVector bfield = pSetup_->inTesla(phase2det->surface().position());
	  LogDebug("PixelDigitizer") << "B-field(T) at " << phase2det->surface().position() << "(cm): " 
	         		     << pSetup_->inTesla(phase2det->surface().position());
	  if (algomap_.find(algotype) != algomap_.end()) 
	    algomap_[algotype]->accumulateSimHits(it, itEnd, globalSimHitIndex, tofBin, phase2det, bfield);
	  else
	    edm::LogInfo("Phase2TrackerDigitizer") << "Unsupported algorithm: " << algotype;
	}
      }
    }
  }
  
  void
  Phase2TrackerDigitizer::initializeEvent(edm::Event const& e, edm::EventSetup const& iSetup) {
    
    // Must initialize all the algorithms
    for (auto it = algomap_.begin(); it != algomap_.end(); ++it) {
      if (first_) it->second->init(iSetup); 
      it->second->initializeEvent(); 
    }
    first_ = false;
    // Make sure that the first crossing processed starts indexing the sim hits from zero.
    // This variable is used so that the sim hits from all crossing frames have sequential
    // indices used to create the digi-sim link (if configured to do so) rather than starting
    // from zero for each crossing.
    crossingSimHitIndexOffset_.clear();
  
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryType_, pDD_);
    iSetup.get<IdealMagneticFieldRecord>().get(pSetup_);
    iSetup.get<IdealGeometryRecord>().get(tTopoHand);
    
    // FIX THIS! We only need to clear and (re)fill this map when the geometry type IOV changes.  Use ESWatcher to determine this.
    if (true) { // Replace with ESWatcher 
      detectorUnits_.clear();
      for (auto iu = pDD_->detUnits().begin(); iu != pDD_->detUnits().end(); ++iu) {
	unsigned int detId_raw = (*iu)->geographicalId().rawId();
	DetId detId = DetId(detId_raw);
	if (DetId(detId).det() == DetId::Detector::Tracker) {
	  Phase2TrackerGeomDetUnit* pixdet = dynamic_cast<Phase2TrackerGeomDetUnit*>(*iu);
	  assert(pixdet);
	  detectorUnits_.insert(std::make_pair(detId_raw, pixdet));
	}
      }
    }
  }  
  void
  Phase2TrackerDigitizer::accumulate(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
    accumulate_local<edm::Event>(iEvent, iSetup);
  }

  void
  Phase2TrackerDigitizer::accumulate(PileUpEventPrincipal const& iEvent, edm::EventSetup const& iSetup) {
    accumulate_local<PileUpEventPrincipal>(iEvent, iSetup);
  }

  template <class T>
  void Phase2TrackerDigitizer::accumulate_local(T const& iEvent, edm::EventSetup const& iSetup) {
    for (auto it = trackerContainers_.begin(); it != trackerContainers_.end(); ++it) {
      edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer_, *it);

      iEvent.getByLabel(tag, simHits);
      unsigned int tofBin = PixelDigiSimLink::LowTof;
      if ((*it).find(std::string("HighTof")) != std::string::npos) tofBin = PixelDigiSimLink::HighTof;
      accumulatePixelHits(simHits, crossingSimHitIndexOffset_[tag.encode()], tofBin);
      // Now that the hits have been processed, I'll add the amount of hits in this crossing on to
      // the global counter. Next time accumulateStripHits() is called it will count the sim hits
      // as though they were on the end of this collection.
      // Note that this is only used for creating digi-sim links (if configured to do so).
      if (simHits.isValid()) crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
     }
  }
  void Phase2TrackerDigitizer::finalizeEvent(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    addPixelCollection(iEvent, iSetup);
    addOuterTrackerCollection(iEvent, iSetup);
  }
  void Phase2TrackerDigitizer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {
  }
  std::string Phase2TrackerDigitizer::getAlgoType(unsigned int detId_raw) {
    DetId detId(detId_raw); 

    std::string algotype = "";
    switch(pDD_->getDetectorType(detId)){
    
    case TrackerGeometry::ModuleType::Ph1PXB:
      algotype = InnerPixel;
      break;
    case TrackerGeometry::ModuleType::Ph1PXF:
      algotype = InnerPixel;
      break;
    case TrackerGeometry::ModuleType::Ph2PSP:
      algotype = PixelinPS;
      break;
    case TrackerGeometry::ModuleType::Ph2PSS:
      algotype = StripinPS;
      break;
    case TrackerGeometry::ModuleType::Ph2SS:
      algotype = TwoStrip;
      break;
    default:
      edm::LogError("Phase2TrackerDigitizer")<<"ERROR - Wrong Detector Type, No Algorithm available "<< pDD_->getDetectorType(detId_raw);
    }
 
    return algotype;
  }
  void Phase2TrackerDigitizer::addPixelCollection(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const TrackerTopology* tTopo = tTopoHand.product();
    std::vector<edm::DetSet<PixelDigi> > digiVector;
    std::vector<edm::DetSet<PixelDigiSimLink> > digiLinkVector;
    for (auto iu = pDD_->detUnits().begin(); iu != pDD_->detUnits().end(); ++iu) {
      DetId detId_raw = DetId((*iu)->geographicalId().rawId());
      const std::string algotype = getAlgoType(detId_raw);

      if (algomap_.find(algotype) == algomap_.end() || algotype != Phase2TrackerDigitizer::InnerPixel) continue;
      std::map<int, DigitizerUtility::DigiSimInfo> digi_map;
      algomap_[algotype]->digitize(dynamic_cast<Phase2TrackerGeomDetUnit*>((*iu)),
                                   digi_map,tTopo);
      edm::DetSet<PixelDigi> collector((*iu)->geographicalId().rawId());
      edm::DetSet<PixelDigiSimLink> linkcollector((*iu)->geographicalId().rawId());
      for (auto it = digi_map.begin(); it != digi_map.end(); ++it) {
	DigitizerUtility::DigiSimInfo info = it->second;  
	std::pair<int,int> ip = PixelDigi::channelToPixel(it->first);
	collector.data.emplace_back(ip.first, ip.second, info.sig_tot);
        for (auto jt = info.track_map.begin(); jt != info.track_map.end(); jt++) {
	  linkcollector.data.emplace_back(it->first, jt->first, info.hit_counter, info.tof_bin, info.event_id, jt->second);
	}
      }  	
      if (collector.data.size() > 0) digiVector.push_back(std::move(collector));	  
      if (linkcollector.data.size() > 0) digiLinkVector.push_back(std::move(linkcollector));
    } 
    
    // Step C: create collection with the cache vector of DetSet 
    std::auto_ptr<edm::DetSetVector<PixelDigi> > 
      output(new edm::DetSetVector<PixelDigi>(digiVector));
    std::auto_ptr<edm::DetSetVector<PixelDigiSimLink> > 
      outputlink(new edm::DetSetVector<PixelDigiSimLink>(digiLinkVector));
    
    // Step D: write output to file 
    iEvent.put(output, "Pixel");
    iEvent.put(outputlink, "Pixel");
  }
  void Phase2TrackerDigitizer::addOuterTrackerCollection(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    const TrackerTopology* tTopo = tTopoHand.product();
    std::vector<edm::DetSet<Phase2TrackerDigi> > digiVector;
    std::vector<edm::DetSet<PixelDigiSimLink> > digiLinkVector;
    for (auto iu = pDD_->detUnits().begin(); iu != pDD_->detUnits().end(); ++iu) {
      DetId detId_raw = DetId((*iu)->geographicalId().rawId());
      const std::string algotype = getAlgoType(detId_raw);

      if (algomap_.find(algotype) == algomap_.end() || algotype == Phase2TrackerDigitizer::InnerPixel) continue;

      std::map<int, DigitizerUtility::DigiSimInfo> digi_map;
      algomap_[algotype]->digitize(dynamic_cast<Phase2TrackerGeomDetUnit*>((*iu)),
				   digi_map, tTopo);
      edm::DetSet<Phase2TrackerDigi> collector((*iu)->geographicalId().rawId());
      edm::DetSet<PixelDigiSimLink> linkcollector((*iu)->geographicalId().rawId());

      for (auto it = digi_map.begin(); it != digi_map.end(); ++it) {
	DigitizerUtility::DigiSimInfo info = it->second;  
	std::pair<int,int> ip = Phase2TrackerDigi::channelToPixel(it->first);
	collector.data.emplace_back(ip.first, ip.second);
        for (auto jt = info.track_map.begin(); jt != info.track_map.end(); jt++) {
	  linkcollector.data.emplace_back(it->first, jt->first, info.hit_counter, info.tof_bin, info.event_id, jt->second);
	}
      }  	
	
      if (collector.data.size() > 0) digiVector.push_back(std::move(collector));	  
      if (linkcollector.data.size() > 0) digiLinkVector.push_back(std::move(linkcollector));
    } 
    
    // Step C: create collection with the cache vector of DetSet 
    std::auto_ptr<edm::DetSetVector<Phase2TrackerDigi> > 
      output(new edm::DetSetVector<Phase2TrackerDigi>(digiVector));
    std::auto_ptr<edm::DetSetVector<PixelDigiSimLink> > 
      outputlink(new edm::DetSetVector<PixelDigiSimLink>(digiLinkVector));
    
    // Step D: write output to file 
    iEvent.put(output, "Tracker");
    iEvent.put(outputlink, "Tracker");
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"

using cms::Phase2TrackerDigitizer;
DEFINE_DIGI_ACCUMULATOR(Phase2TrackerDigitizer);
