#\!/bin/bash
/gridsoft/ipnls/parallel/bin/parallel -j 24 < /gridgroup/cms/viret/SLHC/CMSSW_10_0_0_pre3/src/L1Trigger/TrackFindingAM/test/batch/run_PRMERGE_7_100_test_ele.sh
/gridsoft/ipnls/parallel/bin/parallel -j 24 < /gridgroup/cms/viret/SLHC/CMSSW_10_0_0_pre3/src/L1Trigger/TrackFindingAM/test/batch/run_FMERGE_7_100_test_ele.sh
