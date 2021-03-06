#include "PixelTrackProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include <vector>

using namespace pixeltrackfitting;
using edm::ParameterSet;

PixelTrackProducer::PixelTrackProducer(const ParameterSet& cfg)
  : theReconstruction(cfg)
{
  edm::LogInfo("PixelTrackProducer")<<" construction...";
  produces<reco::TrackCollection>();
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
}

PixelTrackProducer::~PixelTrackProducer() { }

void PixelTrackProducer::endRun(const edm::Run &run, const edm::EventSetup& es)
{ 
  theReconstruction.halt();
}

void PixelTrackProducer::beginRun(const edm::Run &run, const edm::EventSetup& es)
{
  theReconstruction.init(es);
}

void PixelTrackProducer::produce(edm::Event& ev, const edm::EventSetup& es)
{
  LogDebug("PixelTrackProducer, produce")<<"event# :"<<ev.id();

  TracksWithTTRHs tracks;
  theReconstruction.run(tracks,ev,es);

  edm::ESHandle<TrackerTopology> httopo;
  es.get<IdealGeometryRecord>().get(httopo);

  // store tracks
  store(ev, tracks, *httopo);
}

void PixelTrackProducer::store(edm::Event& ev, const TracksWithTTRHs& tracksWithHits, const TrackerTopology& ttopo)
{
  std::auto_ptr<reco::TrackCollection> tracks(new reco::TrackCollection());
  std::auto_ptr<TrackingRecHitCollection> recHits(new TrackingRecHitCollection());
  std::auto_ptr<reco::TrackExtraCollection> trackExtras(new reco::TrackExtraCollection());

  int cc = 0, nTracks = tracksWithHits.size();

  for (int i = 0; i < nTracks; i++)
  {
    reco::Track* track =  tracksWithHits.at(i).first;
    const SeedingHitSet& hits = tracksWithHits.at(i).second;

    for (unsigned int k = 0; k < hits.size(); k++)
    {
      TrackingRecHit *hit = hits[k]->hit()->clone();

      track->setHitPattern(*hit, k, ttopo);
      recHits->push_back(hit);
    }
    tracks->push_back(*track);
    delete track;

  }

  LogDebug("TrackProducer") << "put the collection of TrackingRecHit in the event" << "\n";
  edm::OrphanHandle <TrackingRecHitCollection> ohRH = ev.put( recHits );


  for (int k = 0; k < nTracks; k++)
  {
    reco::TrackExtra* theTrackExtra = new reco::TrackExtra();

    //fill the TrackExtra with TrackingRecHitRef
    unsigned int nHits = tracks->at(k).numberOfValidHits();
    for(unsigned int i = 0; i < nHits; ++i) {
      theTrackExtra->add(TrackingRecHitRef(ohRH,cc));
      cc++;
    }

    trackExtras->push_back(*theTrackExtra);
    delete theTrackExtra;
  }

  LogDebug("TrackProducer") << "put the collection of TrackExtra in the event" << "\n";
  edm::OrphanHandle<reco::TrackExtraCollection> ohTE = ev.put(trackExtras);

  for (int k = 0; k < nTracks; k++)
  {
    const reco::TrackExtraRef theTrackExtraRef(ohTE,k);
    (tracks->at(k)).setExtra(theTrackExtraRef);
  }

  ev.put(tracks);

}
