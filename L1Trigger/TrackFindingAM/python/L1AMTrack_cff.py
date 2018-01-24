import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.MagneticField_cff import *
from L1Trigger.TrackFindingAM.L1AMTrack_cfi import *
from SimTracker.TrackTriggerAssociation.TTStubAssociation_cfi import * 
from SimTracker.TrackTriggerAssociation.TTClusterAssociation_cfi import *
from SimTracker.TrackTriggerAssociation.TTTrackAssociation_cfi import *
import FWCore.ParameterSet.Config as cms


############################################
# STEP 1: AM-based pattern recognition 
############################################

# The simple sequence, creates only a pattern container
# used in principle only for multi-bank PR
#
# Indeed in this case we create the filtered stub container after merging all
# the pattern container (merging the filtered stub container is not possible
# due to persistency loss) 

TTPatternsFromStubs   = cms.Sequence(TTPatternsFromStub)


# The complete sequence, creates the pattern container and the 
# container of filtered stubs/clusters, with corresponding
# associators containers

TTPatternsFromStubswStubs   = cms.Sequence(TTPatternsFromStub*MergePROutput*TTStubAssociatorFromPixelDigis)

############################################
# STEP 2: Combination builder 
############################################

# The simple sequence, creates only a track container
# used in principle only for debugging purposes
#


TTTCsFromPatterns   = cms.Sequence(PatternFilter*TTTCsFromPattern)
TTACBsFromPatterns  = cms.Sequence(PatternFilter*TTACBsFromPattern)
TTSCBsFromPatterns  = cms.Sequence(PatternFilter*TTSCBsFromPattern)

# The sequence. Note that we call the Merge plugins because the filtered containers are created
# here. We just merge one branch...

TTTCsFromPatternswStubs   = cms.Sequence(TTTCsFromPatterns*MergeTCOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)
TTACBsFromPatternswStubs  = cms.Sequence(TTACBsFromPatterns*MergeACBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)
TTSCBsFromPatternswStubs  = cms.Sequence(TTSCBsFromPatterns*MergeSCBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)

# The full combination plan
TTCBsFromPatternswStubs   = cms.Sequence(TTTCsFromPattern*TTACBsFromPattern*TTSCBsFromPattern*MergeTCOutput*MergeACBOutput*MergeSCBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)


############################################
# STEP 3: PCA fit 
############################################

# The simple sequence, creates only a track container
# used in principle only for debugging purposes
#

TTTracksFromTCs   = cms.Sequence(TTTracksTAMUFromTC)
TTTracksFromACBs  = cms.Sequence(TTTracksTAMUFromACB)
TTTracksFromSCBs  = cms.Sequence(TTTracksTAMUFromSCB)



# The sequence. Note that we call the Merge plugins because the filtered containers are created
# here. We just merge one branch...

#TTTracksFromTCswStubs   = cms.Sequence(TTTracksINFNFromTC*MergeFITOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)
TTTracksFromTCswStubs    = cms.Sequence(TTTracksTAMUFromTC*MergeFITOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)
TTTracksFromACBswStubs   = cms.Sequence(TTTracksTAMUFromACB*MergeFITACBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)
TTTracksFromSCBswStubs   = cms.Sequence(TTTracksTAMUFromSCB*MergeFITSCBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)

# The full combination plan
TTTracksFromCBswStubs   = cms.Sequence(TTTracksTAMUFromTC*TTTracksTAMUFromACB*TTTracksTAMUFromSCB*MergeFITOutput*MergeFITACBOutput*MergeFITSCBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)



# Duplicate removal
#DuplicateRemovalTAMUs  = cms.Sequence(DuplicateRemovalTAMU)

DuplicateRemovalswStubs = cms.Sequence(DuplicateRemovalFromTCB*DuplicateRemovalFromACB*DuplicateRemovalFromSCB*MergeDRTCBOutput*MergeDRACBOutput*MergeDRSCBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)

DuplicateRemovalsx2wStubs = cms.Sequence(DuplicateRemovalFromTCB0*DuplicateRemovalFromACB0*DuplicateRemovalFromSCB0*DuplicateRemovalFromTCB*DuplicateRemovalFromACB*DuplicateRemovalFromSCB*MergeDRTCBOutput*MergeDRACBOutput*MergeDRSCBOutput*TTStubAssociatorFromPixelDigis*TTTrackAssociatorFromPixelDigis)


############################################
# STEP 4: MERGE PR outputs
############################################

# This sequence is used mainly the multi-bank merging process, please note that the filtered cluster container is
# not associated due to the lack of simPixelDigis in official samples

TTStubAssociatorFromPixelDigis.TTStubs        = cms.VInputTag( cms.InputTag("MergePROutput", "StubInPattern"))
TTStubAssociatorFromPixelDigis.TTClusterTruth = cms.VInputTag( cms.InputTag("TTClusterAssociatorFromPixelDigis","ClusterAccepted"))

MergePROutputs  = cms.Sequence(MergePROutput*TTStubAssociatorFromPixelDigis)
