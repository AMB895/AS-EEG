#!/usr/bin/env bash
#
# sync to supercomputer
#
# 20250228WF - init
# rsync EEG data for correct AS and error AS trials (prep and stim onset)
# rsync -ihvaz /Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/*cor.{set,fdt} abeatty@bridges2.psc.edu:~/antisaccade_eeg/preprocessed_data/prep/

rsync -ihvaz /Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize_stimonset/*cor.{set,fdt} abeatty@bridges2.psc.edu:~/antisaccade_eeg/preprocessed_data/stim/

# rsync EEG data for VGS trials
# rsync -ihvaz /Volumes/Hera/Abby/preprocessed_data/vgs/AfterWhole/epochclean_homogenize/*kept.{set,fdt} abeatty@bridges2.psc.edu:~/antisaccade_eeg/preprocessed_data/prep/

rsync -ihvaz /Volumes/Hera/Abby/preprocessed_data/vgs/AfterWhole/epochclean_homogenize_stimonset/*kept.{set,fdt} abeatty@bridges2.psc.edu:~/antisaccade_eeg/preprocessed_data/stim/

# rsync calc_ersp_timefreq.m to functions
# rsync -ihvaz /Volumes/Hera/Abby/AS_EEG/erspFunctions/calc_ersp_timefreq.m abeatty@bridges2.psc.edu:~/antisaccade_eeg/functions/

# rsync merge7t_eeg.csv and template channels to PSC
# rsync -ihvaz /Volumes/Hera/Abby/AS_EEG/merge7t_eeg.csv /Volumes/Hera/Abby/AS_EEG/templatechannellabels.mat abeatty@bridges2.psc.edu:~/antisaccade_eeg/
