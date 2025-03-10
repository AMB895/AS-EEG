#!/usr/bin/env bash
#
# sync to supercomputer
#
# 20250228WF - init
#
rsync -ihvaz -n abeatty@bridges2.psc.edu:~/antisaccade_eeg/[p]*.mat /Volumes/Hera/Abby/AS_EEG/TFCE/ 
