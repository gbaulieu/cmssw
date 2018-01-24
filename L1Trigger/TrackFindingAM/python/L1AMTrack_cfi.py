import FWCore.ParameterSet.Config as cms

# AM-based pattern recognition default sequence
TTPatternsFromStub = cms.EDProducer("TrackFindingAMProducer",
   TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTPatternName      = cms.string("AML1Patterns"),
   inputBankFile      = cms.string('/afs/cern.ch/work/s/sviret/testarea/PatternBanks/BE_5D/Eta7_Phi8/ss32_cov40/612_SLHC6_MUBANK_lowmidhig_sec37_ss32_cov40.pbk'),
   threshold          = cms.int32(5),
   nbMissingHits      = cms.int32(-1),
   flatBarrel         = cms.bool(False)
)

## Bend based filtering stage

PatternFilter = cms.EDProducer("RoadCleaningProducer",
                               TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                               TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                               TTTrackName        = cms.string("AML1Patterns"),
                               ptRoadFile         = cms.string("/afs/cern.ch/work/s/sviret/testarea/CMSSW_9_1_0_pre1/src/L1Trigger/TrackFindingAM/data/pt.root"),
                               applyCuts          = cms.bool(False),
                               flatBarrel         = cms.bool(False)
                            )

## Combination building default sequence


TTTCsFromPattern = cms.EDProducer("TrackFitTCProducer",
                                  TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                  TTInputPatterns    = cms.InputTag("PatternFilter", "AML1Patterns"),
#                                  TTInputPatterns    = cms.InputTag("MergePROutput", "AML1Patterns"),
                                  TTTrackName        = cms.string("AML1TCfs"),
                                  TTTrackBinaryName  = cms.string("AML1TCs"),
                                  maxStubsperRoad    = cms.int32(-1),
                                  maxSeedsperTC      = cms.int32(-1),
                                  maxRoadsperTower   = cms.int32(-1),
                                  maxTCsperTower     = cms.int32(-1),
                                  flatBarrel         = cms.bool(False)
                                  )
                     

#TTSCBsFromPattern = ( cms.EDProducer("TrackFitCBProducer",
TTSCBsFromPattern = ( cms.EDProducer("TrackFitPDDSProducer",
                                     TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                     TTInputPatterns    = cms.InputTag("PatternFilter", "AML1Patterns"),
                                     TTTrackName        = cms.string("AML1SCBs"),
                                     maxRoadsperTower   = cms.int32(-1),
                                     maxFitsperTower    = cms.int32(-1),
                                     AdvancedCombinationBuilder = cms.bool(False)
                                     )
                     )

#TTACBsFromPattern = ( cms.EDProducer("TrackFitCBProducer",
TTACBsFromPattern = ( cms.EDProducer("TrackFitPDDSProducer",
                                     TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                     TTInputPatterns    = cms.InputTag("PatternFilter", "AML1Patterns"),
                                     TTTrackName        = cms.string("AML1ACBs"),
                                     maxRoadsperTower   = cms.int32(-1),
                                     maxFitsperTower    = cms.int32(-1),
			             AdvancedCombinationBuilder = cms.bool(True)
                                     )
                      )


## Track fit default sequence

TTTracksINFNFromTC = ( cms.EDProducer("TrackFitPCAProducer",
                                  TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                  TTInputPatterns    = cms.InputTag("MergeTCOutput", "AML1TCfs"),
                                  TTTrackName        = cms.string("AML1Tracks"),
                                  verboseLevel       = cms.untracked.int32(1),
                                  fitPerTriggerTower = cms.untracked.bool(False),
                                  removeDuplicates   = cms.untracked.int32(1)
                                  )
                   )

TTTracksTAMUFromTC = ( cms.EDProducer("AMTrackProducer",
                                  TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                  TTInputPatterns    = cms.InputTag("MergeTCOutput", "AML1TCfs"),
                                  TTTrackName        = cms.string("AMTCL1Tracks"),
                                  CutOnPrincipals    = cms.bool(True)
                                  )
                   )

TTTracksTAMUFromSCB = ( cms.EDProducer("AMTrackProducer",
                                       TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                       TTInputPatterns    = cms.InputTag("MergeSCBOutput", "AML1SCBs"),
                                       TTTrackName        = cms.string("AMSCBL1Tracks"),
                                       CutOnPrincipals    = cms.bool(True)
                                       )
                        )

TTTracksTAMUFromACB = ( cms.EDProducer("AMTrackProducer",
                                       TTInputStubs       = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
                                       TTInputPatterns    = cms.InputTag("MergeACBOutput", "AML1ACBs"),
                                       TTTrackName        = cms.string("AMACBL1Tracks"),
                                       CutOnPrincipals    = cms.bool(True)
                                       )
                        )


#Duplicate removal

DuplicateRemovalFromTCB0 = ( cms.EDProducer("TrackFitDRProducer",
                                           TTTrackName        = cms.InputTag("MergeFITOutput", "AMTCL1Tracks"),
                                           CleanedTTTrackName = cms.string("AMTCL1Tracks"),
                                           ParameterBased     = cms.bool(False),
                                           MaxCommonStubs     = cms.uint32(2)
                                           )
                            )

DuplicateRemovalFromSCB0 = ( cms.EDProducer("TrackFitDRProducer",
                                           TTTrackName        = cms.InputTag("MergeFITSCBOutput", "AMSCBL1Tracks"),
                                           CleanedTTTrackName = cms.string("AMSCBL1Tracks"),
                                           ParameterBased     = cms.bool(False),
                                           MaxCommonStubs     = cms.uint32(2)
                                           )
                            )

DuplicateRemovalFromACB0 = ( cms.EDProducer("TrackFitDRProducer",
                                           TTTrackName        = cms.InputTag("MergeFITACBOutput", "AMACBL1Tracks"),
                                           CleanedTTTrackName = cms.string("AMACBL1Tracks"),
                                           ParameterBased     = cms.bool(False),
                                           MaxCommonStubs     = cms.uint32(2)
                                           )
                            )

DuplicateRemovalFromTCB = ( cms.EDProducer("TrackFitDRProducer",
                                           TTTrackName        = cms.InputTag("DuplicateRemovalFromTCB0", "AMTCL1Tracks"),
                                           CleanedTTTrackName = cms.string("AMTCL1TracksDR"),
                                           ParameterBased     = cms.bool(False),
                                           MaxCommonStubs     = cms.uint32(2)
                                           )
                            )

DuplicateRemovalFromSCB = ( cms.EDProducer("TrackFitDRProducer",
                                           TTTrackName        = cms.InputTag("DuplicateRemovalFromSCB0", "AMSCBL1Tracks"),
                                           CleanedTTTrackName = cms.string("AMSCBL1TracksDR"),
                                           ParameterBased     = cms.bool(False),
                                           MaxCommonStubs     = cms.uint32(2)
                                           )
                            )

DuplicateRemovalFromACB = ( cms.EDProducer("TrackFitDRProducer",
                                           TTTrackName        = cms.InputTag("DuplicateRemovalFromACB0", "AMACBL1Tracks"),
                                           CleanedTTTrackName = cms.string("AMACBL1TracksDR"),
                                           ParameterBased     = cms.bool(False),
                                           MaxCommonStubs     = cms.uint32(2)
                                           )
                            )


# AM output merging sequence
MergePROutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTPatternsFromStub", "AML1Patterns")),                               
   TTFiltClustersName  = cms.string("ClusInPattern"),
   TTFiltStubsName     = cms.string("StubInPattern"),
   TTPatternsName      = cms.string("AML1Patterns")                         
)

MergeTCOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTCsFromPattern", "AML1TCfs")),                               
   TTFiltClustersName  = cms.string("ClusInTC"),
   TTFiltStubsName     = cms.string("StubInTC"),
   TTPatternsName      = cms.string("AML1TCfs")                         
)

MergeACBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTACBsFromPattern", "AML1ACBs")),                               
   TTFiltClustersName  = cms.string("ClusInACB"),
   TTFiltStubsName     = cms.string("StubInACB"),
   TTPatternsName      = cms.string("AML1ACBs")                         
)

MergeSCBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTSCBsFromPattern", "AML1SCBs")),                               
   TTFiltClustersName  = cms.string("ClusInSCB"),
   TTFiltStubsName     = cms.string("StubInSCB"),
   TTPatternsName      = cms.string("AML1SCBs")                         
)


MergeTCOutputb = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTCsFromPattern", "AML1TCfs")),                               
   TTFiltClustersName  = cms.string("ClusInTC"),
   TTFiltStubsName     = cms.string("StubInTC"),
   TTPatternsName      = cms.string("AML1TCfs")                         
)

MergeFITOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksTAMUFromTC", "AMTCL1Tracks")),                               
   TTFiltClustersName  = cms.string("ClusInTCTrack"),
   TTFiltStubsName     = cms.string("StubInTCTrack"),
   TTPatternsName      = cms.string("AMTCL1Tracks")                         
)

MergeFITACBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksTAMUFromACB", "AMACBL1Tracks")),                               
   TTFiltClustersName  = cms.string("ClusInACBTrack"),
   TTFiltStubsName     = cms.string("StubInACBTrack"),
   TTPatternsName      = cms.string("AMACBL1Tracks")                         
)

MergeFITSCBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("TTTracksTAMUFromSCB", "AMSCBL1Tracks")),                               
   TTFiltClustersName  = cms.string("ClusInSCBTrack"),
   TTFiltStubsName     = cms.string("StubInSCBTrack"),
   TTPatternsName      = cms.string("AMSCBL1Tracks")                         
)

MergeDRTCBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("DuplicateRemovalFromTCB", "AMTCL1TracksDR")),
   TTFiltClustersName  = cms.string("ClusInTCBTrackDR"),
   TTFiltStubsName     = cms.string("StubInTCBTrackDR"),
   TTPatternsName      = cms.string("AMTCL1TracksDR")
)

MergeDRACBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("DuplicateRemovalFromACB", "AMACBL1TracksDR")),
   TTFiltClustersName  = cms.string("ClusInACBTrackDR"),
   TTFiltStubsName     = cms.string("StubInACBTrackDR"),
   TTPatternsName      = cms.string("AMACBL1TracksDR")
)

MergeDRSCBOutput = cms.EDProducer("AMOutputMerger",
   TTInputClusters     = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterAccepted"),
   TTInputStubs        = cms.InputTag("TTStubsFromPhase2TrackerDigis", "StubAccepted"),
   TTInputPatterns     = cms.VInputTag(cms.InputTag("DuplicateRemovalFromSCB", "AMSCBL1TracksDR")),
   TTFiltClustersName  = cms.string("ClusInSCBTrackDR"),
   TTFiltStubsName     = cms.string("StubInSCBTrackDR"),
   TTPatternsName      = cms.string("AMSCBL1TracksDR")
)
