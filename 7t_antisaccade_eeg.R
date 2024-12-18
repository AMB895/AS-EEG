merge7t <- read.csv('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv')
antiBehav <- merge7t[c("lunaid","eeg.date","visitno","eeg.age","antiET.Err","antiET.Cor","antiET.ErrCor","antiET.Dropped","antiET.total")]

# calculating error rate for each subject at each time point
error_rate <- antiBehav$antiET.Err/(antiBehav$antiET.Err+antiBehav$antiET.Cor+antiBehav$antiET.ErrCor)
error_rate_bea <- antiBehav$antiET.Err/(antiBehav$antiET.Cor+antiBehav$antiET.ErrCor)
subj_error_rate <- data.frame(antiBehav$lunaid,antiBehav$visitno,antiBehav$eeg.age,error_rate,error_rate_bea)
