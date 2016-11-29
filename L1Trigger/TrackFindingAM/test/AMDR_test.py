#########################
#
# Configuration file for L1 PCA fit
# using a file with AMTC output 
#
# This script works on any official production sample
# (assuming that this sample contains a container of TTStubs,
# a container of TTClusters, and a container of TrackingParticles)
#
# And of course, a container of TCs.... (TTTracks) 
#
#
# Author: S.Viret (viret@in2p3.fr)
# Date        : 04/03/2016
#
# Script tested with release CMSSW_6_2_0_SLHC27
#
#########################

import FWCore.ParameterSet.Config as cms

process = cms.Process('AMDR')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10D_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('L1Trigger.TrackFindingAM.L1AMTrack_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
#
# You can use as input file the result of the script AMTC_test.py of part 6.1.2 of the tutorial
#
# Any other EDM file containing TCs and produced with CMSSW 620_SLHC27 should also work
#

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('file:/fdata/hepx/store/user/demattia/SebsFiles_Tracks/TTbar_PU140_ACB_FIT_output_38.root'),
                            fileNames = cms.untracked.vstring('file:AMFIT_output.root'),
                            duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
)

# The name of the stub container over which the association is done, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

process.TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergeDROutput", "StubInTrack"))
process.TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))
process.TTTrackAssociatorFromPixelDigis.TTTracks      = cms.VInputTag( cms.InputTag("MergeDROutput", "AML1Tracks"))

# process.TTTracksTAMUFromCB.ConstantsDir               = cms.FileInPath("L1Trigger/TrackFindingAM/data/PreEstimate_Transverse/matrixVD_2016.txt")

# Additional output definition
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('TTbar_PU140_ACB_DR_output_38_parameterBased.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

# For the moment need to explicitely keep the following containers
# (not yet in the customizing scripts)

# Keep the PR output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMPR')
# process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMTC')
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMCB')
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMFIT')

# Keep the FIT output
process.RAWSIMoutput.outputCommands.append('keep  *_*_*_AMDR')
process.RAWSIMoutput.outputCommands.append('drop *_TTTracks*FromCB_*_*')
process.RAWSIMoutput.outputCommands.append('keep  *_*_MergedTrackTruth_*')



process.RAWSIMoutput.outputCommands.append('keep  *_*_*_*')



# Path and EndPath definitions
process.L1AMDR_step          = cms.Path(process.DuplicateRemovalTAMUswStubs)
process.endjob_step          = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step    = cms.EndPath(process.RAWSIMoutput)

# Uncomment this to use the parameter-based duplicate removal
# process.DuplicateRemovalTAMU.ParameterBased = True

# Use this to change the maximum number of common stubs for the stub-based duplicate removal
# process.DuplicateRemovalTAMU.MaxCommonStubs = 2

process.schedule = cms.Schedule(process.L1AMDR_step,process.endjob_step,process.RAWSIMoutput_step)
# process.schedule = cms.Schedule(process.L1AMCB_step,process.L1AMFIT_step,process.endjob_step,process.RAWSIMoutput_step)


# Automatic addition of the customisation function

from SLHCUpgradeSimulations.Configuration.combinedCustoms import customiseBE5DPixel10D
from SLHCUpgradeSimulations.Configuration.combinedCustoms import customise_ev_BE5DPixel10D

process=customiseBE5DPixel10D(process)
process=customise_ev_BE5DPixel10D(process)

# End of customisation functions
