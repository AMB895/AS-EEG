#!/usr/bin/env bash
#
# sync to supercomputer
#
# 20250228WF - init
#
rsync -ihvaz  PrepPeriodAnalysis/CorrectTrials/[lt]*.mat PrepPeriodAnalysis/*.mat abeatty@bridges2.psc.edu:~/antisaccade_eeg/data/
rsync -ihvaz  /Volumes/Hera/Abby/AS_EEG/limo_tools/ /Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/ /Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/ abeatty@bridges2.psc.edu:~/antisaccade_eeg/tools/
rsync -ihvaz /Volumes/Hera/Abby/AS_EEG/STUDY/logicalchanneighbmatrix.mat abeatty@bridges2.psc.edu:~/antisaccade_eeg/

