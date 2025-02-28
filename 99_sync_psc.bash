#!/usr/bin/env bash
#
# sync to supercomputer
#
# 20250228WF - init
#
rsync -ihvaz  PrepPeriodAnalysis/CorrectTrials/[lt]*.mat PrepPeriodAnalysis/*.mat abeatty@bridges2.psc.edu:~/antisaccade_eeg/data/
