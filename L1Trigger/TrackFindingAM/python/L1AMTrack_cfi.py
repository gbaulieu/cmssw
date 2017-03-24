import FWCore.ParameterSet.Config as cms

# AM-based pattern recognition default sequence
TTPatternsFromStub = cms.EDProducer("TrackFindingAMProducer",
   TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTPatternName      = cms.string("AML1Patterns"),
   inputBankFile      = cms.string('/afs/cern.ch/work/s/sviret/testarea/PatternBanks/BE_5D/Eta7_Phi8/ss32_cov40/612_SLHC6_MUBANK_lowmidhig_sec37_ss32_cov40.pbk'),
   threshold          = cms.int32(5),
   nbMissingHits      = cms.int32(-1)
)

TTTCsFromPattern = ( cms.EDProducer("TrackFitTCProducer",
                                    TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                    TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                    TTTrackName        = cms.string("AML1TCs"),
                                    TTTrackBinaryName  = cms.string("AML1BinTCs")
                                    )
                     )

# AM output merging sequence
MergePROutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTPatternsFromStub", "AML1Patterns")),                               
   TTFiltClustersName  = cms.string("ClusInPattern"),
   TTFiltStubsName     = cms.string("StubInPattern"),
   TTPatternsName      = cms.string("AML1Patterns")                         
)


MergeTCOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTCsFromPattern", "AML1TCs")),                               
   TTFiltClustersName  = cms.string("ClusInTC"),
   TTFiltStubsName     = cms.string("StubInTC"),
   TTPatternsName      = cms.string("AML1TCs")                         
)

MergeTCOutputb = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTStubsFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTCsFromPattern", "AML1BinTCs")),                               
   TTFiltClustersName  = cms.string("ClusInTC"),
   TTFiltStubsName     = cms.string("StubInTC"),
   TTPatternsName      = cms.string("AML1BinTCs")                         
)
