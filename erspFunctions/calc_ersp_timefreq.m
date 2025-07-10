function [ERSP_dB_prep,ERS_allchans,evoked_powerspec_prep,spon_powerspec_prep,times,freqs,preptimes] = calc_ersp_timefreq_WIP(EEG)

% Looping through each channel
numchans = EEG.nbchan;
for currentChan = 1:numchans
    % printing out channel progress
    fprintf('\n-------Channel %d-------\n',currentChan)
    %% Compute ERSP with trial average baseline correction
    % Step 1: Compute spectral estimate at each time point for each trial
    % Inputs to timefreq():
    %   data: 2d array (times,trials)
    %   srate: sampling rate
    %   'cycles': [3 0.5], starts with 3 wavelet cycles and increases linearly
    %   'wletmethod' = 'dftfilt3' (Morlet wavelet/ DEFAULT method)
    %   'tlimits' = epoch min and max times in ms
    %   'timesout' = array of times to obtain spectral decomposition 
    %   'freqs' = frequency limits
    %   'nfreqs'' = number of frequency points
    chanEEGdata = squeeze(EEG.data(currentChan,:,:));
    [tf,freqs,times] = timefreq(chanEEGdata,EEG.srate,'cycles',[3 0.5],'wletmethod','dftfilt3',...
        'tlimits',[EEG.xmin EEG.xmax]*1000,'timesout',-800:7:800,'freqs',[3 70],'nfreqs',150);
    % getting start and end times of baseline period and prep period
    % start and end of baseline
    startbase_idx = find(times>-700,1,'first');
    endbase_idx = find(times>=-400 ,1,'first');
    % start and end of prep
    startprep_idx = find(times>=0,1,'first');
    endprep_idx = find(times>=500,1,'first');
    preptimes = times(startprep_idx:endprep_idx);
    
    % Step 2: Calculate power spectrum of complex spectral estimate for each trial
    powerspec_alltrials = abs(tf).^2;
    
    % From Grandchamp & Delmore, 2011:
    % ERS(f,t) = Average power spectrum for at each time point windows across trials
    % ERSP_%(f,t) = ERS(f,t)/uB(f)
    % ERSP_dB(f,t) = 10*log10(ERSP_%(f,t)) 
    
    % Step 3: Calculate event-related spectrum (ERS)
    ERS = mean(powerspec_alltrials,3); % average power spectrum across trials
    % ^^ also considered total power spectrum for spontaneous power calculations
    
    % Step 4: Define baseline power spectrum (averaged across trials and baseline time points)
    baseline_powerspec = mean(powerspec_alltrials(:,startbase_idx:endbase_idx,:),[2 3]); % averaging across 2nd and 3rd dimensions (baseline times and trials)
    
    % Step 5: Calculate percent signal change from baseline ERSP
    ERSP_percent = ERS./baseline_powerspec;
    
    % Step 6: Convert percent signal change to decibels
    ERSP_dB = 10*log10(ERSP_percent);
    
    % just want to look at prep period ERSP 
    ERSP_dB_prep(currentChan,:,:) = ERSP_dB(:,startprep_idx:endprep_idx);
    
    % save out ERS for each channel
    ERS_allchans(currentChan,:,:) = ERS;
    
    %% Calculate total power
    % total power spectrum = averaging power spectrum across trials (ERS) in frequency domain
    total_powerspec_prep(currentChan,:,:) = ERS(:,startprep_idx:endprep_idx);
    
    %% Calculate evoked power
    % Get evoked power spectrum by averaging trials in time domain then
    % computing power spectrum
    chanEEGdata_avgtrials = squeeze(mean(EEG.data(currentChan,:,:),3))';
    [tf_evoked,freqs,times] = timefreq(chanEEGdata_avgtrials,EEG.srate,'cycles',[3 0.5],'wletmethod','dftfilt3',...
        'tlimits',[EEG.xmin EEG.xmax]*1000,'timesout',-800:7:800,'freqs',[3 70],'nfreqs',150);
    % power across trials averaged in time domain
    powerspec_avgtrials= abs(tf_evoked).^2;
    evoked_powerspec_prep(currentChan,:,:) = powerspec_avgtrials(:,startprep_idx:endprep_idx);
    %% Calculate spontaneous power
    % Total power spectrum - evoked power spectrum
    spon_powerspec_prep(currentChan,:,:) = total_powerspec_prep(currentChan,:,:) - evoked_powerspec_prep(currentChan,:,:);
end
end