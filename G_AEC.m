%% Amplitude envelope correlation using hilbert amplitude envelopes
% amb895 05/07/2026
clear; 

% add paths
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/AS_EEG/')
addpath('/Volumes/Hera/Abby/Resources/slanCM/')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab("nogui");

% load in frequency bands of interest
load('/Volumes/Hera/Abby/AS_EEG/Results/ersp_itc_clusterinfo_avgtrl_20260508.mat')

% just ersp interaction clusters 1 and 5
clust_table = clustertable((clustertable.clusternum == 1 | clustertable.clusternum==5) & ...
    strcmp(clustertable.modelterm, "interaction") & strcmp(clustertable.datatype,"ersp"),:);
betastart = clust_table.minfreq(clust_table.clusternum == 5);
betaend = clust_table.maxfreq(clust_table.clusternum == 5);
gammastart = clust_table.minfreq(clust_table.clusternum == 1);
gammaend = clust_table.maxfreq(clust_table.clusternum == 1);

% load in AS eeg merged 7t table
load('/Volumes/Hera/Abby/AS_EEG/Results/merged7t_ASscore_table.mat')
merged7t_eeg = merged7t_ASscore_table;

% loading template of channel labels
load('/Volumes/Hera/Abby/AS_EEG/templatechannellabels.mat')

% tasks
task = ["corAS","VGS"];

% electrodes
frow = [4 5 6 7 37 38 39 40 41];
nchans = length(frow);

%% Directories
% Main directory of preprocessed data
maindir = hera('Abby/preprocessed_data');

% where to save AEC 
savedir = '/Volumes/Hera/Abby/AS_EEG/Results/';

% stim or prep onset
eventtype = 'prep';

% Directories to AS and VGS preprocessed data
taskdirAS = [maindir '/' 'anti']; 
taskdirVGS = [maindir '/' 'vgs']; 

% Directories to AS and VGS epoched data
% epochdirAS = [taskdirAS,'/AfterWhole/kept_epochclean_prep/'];
% epochdirVGS = [taskdirVGS,'/AfterWhole/kept_epochclean_prep/'];

% % Directories to AS and VGS data (not epoched)
epochdirAS = [taskdirAS,'/AfterWhole/ICAwholeClean_homogenize/'];
epochdirVGS = [taskdirVGS,'/AfterWhole/ICAwholeClean_homogenize/'];

% Epoch data file names
% epochcorAS_filename = ['*_epochs_kept_cor_' eventtype '.set'];
% epochVGS_filename = ['*_epochs_kept_' eventtype '.set'];

%  data file names
epochcorAS_filename = ['*_icapru.set'];
epochVGS_filename = ['*_icapru.set'];

% Add all eeg files to one big structure to loop over
EEGfilenames = {dir([epochdirAS,epochcorAS_filename]) dir([epochdirVGS,epochVGS_filename])};
tasks = ["corAS","VGS"];

% event start and end times (in ms)
eventstart = 0;
eventend = 500;

%% Loop through each task and subject
temptable = merged7t_eeg;
temptable.corAS_AEC = NaN(height(temptable),1);
temptable.corAS_PEC= NaN(height(temptable),1);
temptable.VGS_AEC = NaN(height(temptable),1);
temptable.VGS_PEC = NaN(height(temptable),1);

for t = 1:length(tasks)
    % Defining current task EEG
    taskfilenames = EEGfilenames{t};
    tasktotalEEGs = length(taskfilenames);
    currtask = tasks(t);
    for currentEEG = 1:tasktotalEEGs
        % getting subject's id and eeg date
        filename = [taskfilenames(currentEEG).name];
        filepath = [taskfilenames(currentEEG).folder];
        splitfilename = strsplit(filename,'_');
        lunaid = cell2mat(splitfilename(1));
        eeg_date = cell2mat(splitfilename(2));
        
        % defining subject's merged7t table
        subIdx = find(merged7t_eeg.lunaid==str2double(lunaid) & merged7t_eeg.eeg_date==str2double(eeg_date));
        subTable = merged7t_eeg(subIdx,:);
        
        % skip if subject is not in merged 7t
        if isempty(subTable)
            warning('Subject %s %s is not in merge 7t; skipping\n',lunaid,eeg_date)
            temptable(subIdx,sprintf("%s_AEC",currtask)) = table(NaN);
            temptable(subIdx,sprintf("%s_PEC",currtask)) = table(NaN);            
            continue
        end
        % skip if age is empty in merged 7t table
        if isnan(subTable.eeg_age)
            warning('Subject %s %s does not have age; skipping\n',lunaid,eeg_date)
            temptable(subIdx,sprintf("%s_AEC",currtask)) = table(NaN);
            temptable(subIdx,sprintf("%s_PEC",currtask)) = table(NaN);    
            continue
        end
        % load subject's EEG
        EEG = pop_loadset('filename',filename,'filepath',filepath);
        chanNames = string({EEG.chanlocs.labels}');
        % skipping if channel labels do not match template
        if all(tempChanLabels~=chanNames)
            warning('Check Subject Channel Labels %s %s; skipping\n',lunaid,eeg_date)
            temptable(subIdx,sprintf("%s_AEC",currtask)) = table(NaN);
            temptable(subIdx,sprintf("%s_PEC",currtask)) = table(NaN);    
            continue
        end
        % skip if subject has 0 or 1 epochs (fails if subject has only 1 epoch)
%         numEpochs = length(EEG.epoch);
%         if numEpochs < 1
%             warning('Subject %s %s has %d epochs; skipping\n',lunaid,eeg_date,numEpochs)
%             temptable(subIdx,sprintf("%s_AEC",currtask)) = table(NaN);
%             temptable(subIdx,sprintf("%s_PEC",currtask)) = table(NaN);    
%             continue
%         end
        freqbands = [betastart betaend; gammastart gammaend];
        freqbandsstr = ["beta","gamma"];
        [EEGfilt_epoch, EEGdatafilt_epoch, prepIDX] = BPfilt_epoch(EEG, freqbands, freqbandsstr,currtask);
        
        if isnan(EEGdatafilt_epoch)
            warning('Subject %s %s has 0 or 1 epochs; skipping\n',lunaid,eeg_date)
            temptable(subIdx,sprintf("%s_AEC",currtask)) = table(NaN);
            temptable(subIdx,sprintf("%s_PEC",currtask)) = table(NaN);    
            continue
        end
        % compute amplitude envelope and power envelope
        times = EEGfilt_epoch.beta.times(prepIDX(1):prepIDX(2));
%         beta_amp_env = squeeze(abs(EEGdatafilt_epoch(1,frow,prepIDX(1):prepIDX(2),:)));
%         beta_pow_env = squeeze(abs(EEGdatafilt_epoch(1,frow,prepIDX(1):prepIDX(2),:)).^2);
%         gamma_amp_env =squeeze( abs(EEGdatafilt_epoch(2,frow,prepIDX(1):prepIDX(2),:)));
%         gamma_pow_env = squeeze(abs(EEGdatafilt_epoch(2,frow,prepIDX(1):prepIDX(2),:)).^2);
%         
%         for chan = 1:nchans
%             beta_amp_ch = squeeze(beta_amp_env(chan,:,:));
%             gamma_amp_ch = squeeze(gamma_amp_env(chan,:,:));
%             beta_power_ch = squeeze(beta_pow_env(chan,:,:));
%             gamma_power_ch = squeeze(gamma_pow_env(chan,:,:));
%             ch_AEC(chan) = mean(diag(corr(beta_amp_ch,gamma_amp_ch)));
%             ch_PEC(chan) = mean(diag(corr(beta_power_ch, gamma_power_ch)));
%         end
%         AEC = mean(ch_AEC);
%         PEC = mean(ch_PEC);
        beta_signal = squeeze(mean(EEGdatafilt_epoch(1,frow,prepIDX(1):prepIDX(2),:),2));
        gamma_signal = squeeze(mean(EEGdatafilt_epoch(2,frow,prepIDX(1):prepIDX(2),:),2));
        beta_hilenv = abs(hilbert(beta_signal));
        gamma_hilenv = abs(hilbert(gamma_signal));
        r = corr(beta_hilenv,gamma_hilenv);
        AEC = mean(diag(r));


%         beta_amp_env = squeeze(mean(abs(EEGdatafilt_epoch(1,frow,prepIDX(1):prepIDX(2),:)),[2]));
%         beta_pow_env = squeeze(mean(abs(EEGdatafilt_epoch(1,frow,prepIDX(1):prepIDX(2),:)).^2,[2]));
%         gamma_amp_env = squeeze(mean(abs(EEGdatafilt_epoch(2,frow,prepIDX(1):prepIDX(2),:)),[2]));
%         gamma_pow_env = squeeze(mean(abs(EEGdatafilt_epoch(2,frow,prepIDX(1):prepIDX(2),:)).^2,[2 ]));
%         AEC = mean(diag(corr(beta_amp_env,gamma_amp_env)));
%         PEC = mean(diag(corr(beta_pow_env, gamma_pow_env)));
%         AEC = corr(beta_amp_env,gamma_amp_env);
%         PEC = corr(beta_pow_env, gamma_pow_env);
        % plot beta and gamma
%         figure;
%         subplot(2,1,1)
%         plot(times,beta_amp_env); hold on; plot(times,gamma_amp_env); legend('beta','gamma')
%         
%         subplot(2,1,2)
%         plot(times,beta_pow_env); hold on; plot(times,gamma_pow_env); legend('beta','gamma')
        
        

        % compute tf reconstruction of amplitude and power
%         [power, times, freqs, amp_env] = F_tfpower(EEG,'amplitude',1, 'channels',frow);
%         starttime_idx = find(times>=eventstart,1,'first');
%         endtime_idx = find(times>=eventend,1,'first');
% %         power = squeeze(mean(power(:,starttime_idx:endtime_idx,:,:),[4]));        
% %         amp_env = squeeze(mean(amp_env(:,starttime_idx:endtime_idx,:,:),[ 4] ));
%         power = power(:,starttime_idx:endtime_idx,:,:);        
%         amp_env = amp_env(:,starttime_idx:endtime_idx,:,:);
%         times = times(starttime_idx:endtime_idx);
%         ntrl = size(amp_env,3);
%         betastart_idx = find(freqs>=betastart,1,'first');
%         betaend_idx = find(freqs>=betaend,1,'first');
%         gammastart_idx = find(freqs>=gammastart,1,'first');
%         gammaend_idx = find(freqs>=gammaend,1,'first');
%         betapower = squeeze(mean(power(betastart_idx:betaend_idx,:,:,:),[1 3 4]));
%         betaamp = squeeze(mean(amp_env(betastart_idx:betaend_idx,:,:,:),[1 3 4]));
%         gammapower = squeeze(mean(power(gammastart_idx:gammaend_idx,:,:,:),[1 3 4]));
%         gammaamp = squeeze(mean(amp_env(gammastart_idx:gammaend_idx,:,:,:),[1 3 4]));
        
%         r_a = diag(corr(betaamp,gammaamp));
%         r_p = diag(corr(betapower,gammapower));
        
%         AEC = mean(r_a);
%         PEC = mean(r_p);
        
%         AEC =corr(betaamp',gammaamp');
%         PEC = corr(betapower',gammapower');
%         for trl = 1:ntrl
%             trl_amp = squeeze(amp_env(:,:,trl));
%             trl_power = squeeze(power(:,:,trl));
%             trl_rA = corr(trl_amp');
%             trl_rP = corr(trl_power');
%             trl_AEC = mean(trl_rA(betastart_idx:betaend_idx,gammastart_idx:gammaend_idx,:),'all');
%             trl_PEC = mean(trl_rP(betastart_idx:betaend_idx,gammastart_idx:gammaend_idx,:),'all');
%             AEC(trl) = trl_AEC;
%             PEC(trl) = trl_PEC;
%             clear trl_AEC trl_PEC
%         end
%         
% %         
% %         r_A = corr(amp_env');
% %         r_P = corr(power');
%         AEC = mean(AEC,'all');
%         PEC =  mean(PEC,'all');
        % correlation across trials and time
%         ntimes = length(times);
%         ntrls = size(amp_env,3);
%         nfreqs = length(freqs);
%         
%         A = reshape(amp_env, nfreqs, ntimes*ntrls);
%         r_A = corr(A');
%         P = reshape(power, nfreqs, ntimes*ntrls);
%         r_P = corr(P');
        
        % just beta-gamma correlation        
%         betastart_idx = find(freqs>=betastart,1,'first');
%         betaend_idx = find(freqs>=betaend,1,'first');
%         gammastart_idx = find(freqs>=gammastart,1,'first');
%         gammaend_idx = find(freqs>=gammaend,1,'first');
%         AEC = mean(r_A(betastart_idx:betaend_idx,gammastart_idx:gammaend_idx),'all');
%         PEC = mean(r_P(betastart_idx:betaend_idx,gammastart_idx:gammaend_idx),'all');
                
        % add r_power and r_amp to temp table
        temptable(subIdx,sprintf("%s_AEC",currtask)) = table(AEC);
%         temptable(subIdx,sprintf("%s_PEC",currtask)) = table(PEC);
        
    end
end

%% save table
currdate = datetime();
currdate.Format = 'yyyyMMdd';
AECtable = temptable;
save(fullfile(savedir,['AEC_' char(currdate) '.mat']),'AECtable')
writetable(AECtable,fullfile(savedir,['AEC_' char(currdate) '.csv']))