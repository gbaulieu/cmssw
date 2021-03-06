import FWCore.ParameterSet.Config as cms

def customiseDefault(process):
    if hasattr(process,'pfTrack'):
        process.pfTrack.TrajInEvents = cms.bool(True)
    if hasattr(process,'csc2DRecHits'):
        process.csc2DRecHits.readBadChannels = cms.bool(False)

    if hasattr(process,'validation_step'):
        process.validation_step.remove(process.HLTSusyExoValFastSim)
        process.validation_step.remove(process.hltHiggsValidator)

    return process

def customisePhase2(process):
    process=customiseDefault(process)
    process.famosSimHits.MaterialEffects.PairProduction = cms.bool(False)
    process.famosSimHits.MaterialEffects.Bremsstrahlung = cms.bool(False)
    process.famosSimHits.MaterialEffects.MuonBremsstrahlung = cms.bool(False)
    process.famosSimHits.MaterialEffects.EnergyLoss = cms.bool(False)
    process.famosSimHits.MaterialEffects.MultipleScattering = cms.bool(False)
    # keep NI so to allow thickness to be properly treated in the interaction geometry
    process.famosSimHits.MaterialEffects.NuclearInteraction = cms.bool(True)
    process.KFFittingSmootherWithOutlierRejection.EstimateCut = cms.double(50.0)

    return process

