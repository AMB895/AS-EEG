%% Evoked and Spontaneous CFC
% setting up for correct and vgs trials
clear; close all;
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/preprocessed_data/vgs/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/')
addpath('/Volumes/Hera/Abby/AS_EEG/BPfilteredData/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;
%% Correct AS trials
% Main directory of BP filtered EEG data
maindirfiltered = hera('Abby/AS_EEG/BPfilteredData/');
bpdataname = '*corBPfilteredEEG.mat';
% Getting all filtered data files
filterEEGfileinfo = dir([maindirfiltered,bpdataname]);
totalfilterEEGs = size(filterEEGfileinfo,1);
% Main directory of epoched EEG data
maindirepoch = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindirepoch '/' task]; 
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept_cor.set';
% Getting all epoched data files
epochEEGfileinfo = dir([epochedpath,epochedname]);

%% compute correlation matrix for each subject- correct AS trials
noMatch = [];
oneTrial = [];
if exist('corCorrelationMats.mat','file')
    fprintf('Already computed CFC correlation matricies for correct AS trials\n')
else
    for currentEEG = 1:totalfilterEEGs
        % getting file info of epoched data
        epochEEGfilename = [epochEEGfileinfo(currentEEG).name];
        epochSplitString = split(epochEEGfilename,'_');
        epochSubID = cell2mat(epochSplitString(1));
        epochScanDate = cell2mat(epochSplitString(2));
        % getting file info of filtered data
        filterEEGfilename = [filterEEGfileinfo(currentEEG).name];
        filterSplitString = split(filterEEGfilename,'_');
        filterSubID = cell2mat(filterSplitString(1));
        filterScanDate = cell2mat(filterSplitString(2));

        % make sure the epoch EEG and filtered EEG are from the same person
        if strcmp(epochSubID,filterSubID) && strcmp(epochScanDate,filterScanDate)
            fprintf('Epoch EEG and filtered EEG match subject IDs\n')
            % load epoch EEG to get frames per epoch and total trials
            epochEEG = pop_loadset(epochEEGfilename);
            framesPerEpoch = epochEEG.pnts;
            totalTrials = epochEEG.trials;
            % load filtered EEG for band pass filtered EEGs
            subStruct = load(filterEEGfilename);
            filteredStruct = subStruct.corBPfilterdata;
            thetaEEG = filteredStruct.corBPthetaEEG;
            lowbetaEEG = filteredStruct.corBPlowbetaEEG;
            highbetaEEG = filteredStruct.corBPhighbetaEEG;
            
            % average filtered time-courses across F-row
            thetaEEG_frow = squeeze(mean(thetaEEG([4 5 6 7 37 38 39 40 41],:),1));
            lowbetaEEG_frow = squeeze(mean(lowbetaEEG([4 5 6 7 37 38 39 40 41],:),1));
            highbetaEEG_frow = squeeze(mean(highbetaEEG([4 5 6 7 37 38 39 40 41],:),1));
            
            % reshape filtered EEG
            % original filtered EEG data is 1 x (totalTrials*framesPerEpoch)
            totalFrames = length(thetaEEG_frow);
            testNumTrials = totalFrames/framesPerEpoch;
            
            % check that computed number of trials matches actual number of trials
            if testNumTrials == totalTrials
                fprintf('Computed number of trials and total trials match\n')
            else
                fprintf('Computed number of trials and total trials do not match\n')
                noMatch(end+1,:) = [suID,scanDate];
                continue;
            end
            
             % if subject only has one trial reshape() messes up data
            if totalFrames == framesPerEpoch && totalTrials ==1
                fprintf('Subject %s %s only has 1 trial\n',filterSubID,filterScanDate)
                thetaEEG_frow_trials = thetaEEG_frow;
                lowbetaEEG_frow_trials = lowbetaEEG_frow;
                highbetaEEG_frow_trials = highbetaEEG_frow;
                
                % truncate filtered EEGs to start (O ms) and end (500 ms) of prep period
                time0 = find(epochEEG.times == 0);
                time500 = find(epochEEG.times > 499 & epochEEG.times <501,1);
                thetaEEG_frow_trials = thetaEEG_frow_trials(:,time0:time500);
                lowbetaEEG_frow_trials = lowbetaEEG_frow_trials(:,time0:time500);
                highbetaEEG_frow_trials = highbetaEEG_frow_trials(:,time0:time500);
                oneTrial(end+1,:) = [filterSubID,filterScanDate];
            else
                % reshape filtered EEG data to matrix of trials*framesPerEpoch
                thetaEEG_frow_trials = reshape(thetaEEG_frow,[totalTrials,framesPerEpoch]);
                lowbetaEEG_frow_trials = reshape(lowbetaEEG_frow,[totalTrials,framesPerEpoch]);
                highbetaEEG_frow_trials = reshape(highbetaEEG_frow,[totalTrials,framesPerEpoch]);
                
                % truncate filtered EEGs to start (O ms) and end (500 ms) of prep period
                time0 = find(epochEEG.times==0);
                time500 = find(epochEEG.times==500);
                thetaEEG_frow_trials = thetaEEG_frow_trials(:,time0:time500);
                lowbetaEEG_frow_trials = lowbetaEEG_frow_trials(:,time0:time500);
                highbetaEEG_frow_trials = highbetaEEG_frow_trials(:,time0:time500);
            end
            
            %% Spontaneous CFC- Correct AS Trials
            % sliding window and computing correlation in each window
            framesPerWindow = 15; % samples (100 ms)
            overlap = 10; % samples (66.67 ms)
            for currentTrial = 1:totalTrials
                % Generating time windows for each trial
                windowed_thetaEEG = buffer(thetaEEG_frow_trials(currentTrial,:),framesPerWindow,overlap,'nodelay');
                windowed_lowbetaEEG = buffer(lowbetaEEG_frow_trials(currentTrial,:),framesPerWindow,overlap,'nodelay');
                windowed_highbetaEEG = buffer(highbetaEEG_frow_trials(currentTrial,:),framesPerWindow,overlap,'nodelay');
                % buffer() outputs frames per window x # of windows
                % columns are the time windows, rows are the frames in the window
                numWindows = size(windowed_thetaEEG,2);
                
                % In each window, calculating spontaneous CFC for each trial
                for currentWindow = 1:numWindows
                    currentWindow_thetaEEG = nonzeros(windowed_thetaEEG(:,currentWindow));
                    currentWindow_lowbetaEEG = nonzeros(windowed_lowbetaEEG(:,currentWindow));
                    currentWindow_highbetaEEG = nonzeros(windowed_highbetaEEG(:,currentWindow));
                    mat = [currentWindow_thetaEEG,currentWindow_lowbetaEEG,currentWindow_highbetaEEG];
                    subcorR_spontaneous(currentTrial,currentWindow,:,:) = corrcoef(mat);
                end
            end
            % Average subcorR_spontaneous across trials to get average spontaneous correlation for each subject
            corR_spontaneous(currentEEG,:,:,:) = squeeze(mean(subcorR_spontaneous,1));
            
            %% Evoked CFC- Correct AS Trials
            % average filtered EEGs across trials to get evoked correlation
            thetaEEG_frow_avgtrials = mean(thetaEEG_frow_trials,1);
            lowbetaEEG_frow_avgtrials = mean(lowbetaEEG_frow_trials,1);
            highbetaEEG_frow_avgtrials = mean(highbetaEEG_frow_trials,1);
            windowed_thetaEEG_avgtrials = buffer(thetaEEG_frow_avgtrials,framesPerWindow,overlap,'nodelay');
            windowed_lowbetaEEG_avgtrials = buffer(lowbetaEEG_frow_avgtrials,framesPerWindow,overlap,'nodelay');
            windowed_highbetaEEG_avgtrials = buffer(highbetaEEG_frow_avgtrials,framesPerWindow,overlap,'nodelay');
            % buffer() outputs frames per window x # of windows
            % columns are the time windows, rows are the frames in the window
            numWindows = size(windowed_thetaEEG_avgtrials,2);
            
            % evoked CFC at each time window
            for currentWindow = 1:numWindows
                % nonzeros() b/c buffer zero-pads windows if not evenly divisible
                currentWindow_thetaEEG_avgtrials = nonzeros(windowed_thetaEEG_avgtrials(:,currentWindow));
                currentWindow_lowbetaEEG_avgtrials = nonzeros(windowed_lowbetaEEG_avgtrials(:,currentWindow));
                currentWindow_highbetaEEG_avgtrials = nonzeros(windowed_highbetaEEG_avgtrials(:,currentWindow));
                mat = [currentWindow_thetaEEG_avgtrials,currentWindow_lowbetaEEG_avgtrials,currentWindow_highbetaEEG_avgtrials];
                corR_evoked(currentEEG,currentWindow,:,:) = corrcoef(mat);
            end
        else
            fprintf('Epoch EEG and filtered EEG do not match subject IDs; skipping')
        end
    end
    % save evoked and spontaneous correlation matrix
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corCorrelationMats.mat','corR_evoked','corR_spontaneous')
end

%% VGS trials
% Main directory of BP filtered EEG data
maindirfiltered = hera('Abby/AS_EEG/BPfilteredData/');
bpdataname = '*vgsBPfilteredEEG.mat';
% Getting all filtered data files
filterEEGfileinfo = dir([maindirfiltered,bpdataname]);
totalfilterEEGs = size(filterEEGfileinfo,1);
% Main directory of epoched EEG data
maindirepoch = hera('Abby/preprocessed_data');
task = 'vgs';
taskdirectory = [maindirepoch '/' task]; 
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept.set';
% Getting all epoched data files
epochEEGfileinfo = dir([epochedpath,epochedname]);
totalepochEEGs = size(epochEEGfileinfo,1);

%% compute correlation matrix for each subject- VGS trials
noMatch = [];
if exist('vgsCorrelationMats.mat','file')
    fprintf('Already computed CFC correlation matricies for VGS trials\n')
else
    for currentEEG = 1:totalfilterEEGs
        % getting file info of epoched data
        epochEEGfilename = [epochEEGfileinfo(currentEEG).name];
        epochSplitString = split(epochEEGfilename,'_');
        epochSubID = cell2mat(epochSplitString(1));
        epochScanDate = cell2mat(epochSplitString(2));
        % getting file info of filtered data
        filterEEGfilename = [filterEEGfileinfo(currentEEG).name];
        filterSplitString = split(filterEEGfilename,'_');
        filterSubID = cell2mat(filterSplitString(1));
        filterScanDate = cell2mat(filterSplitString(2));

        % make sure the epoch EEG and filtered EEG are from the same person
        if strcmp(epochSubID,filterSubID) && strcmp(epochScanDate,filterScanDate)
            fprintf('Epoch EEG and filtered EEG match subject IDs\n')
            % load epoch EEG for frames per epoch and total trials
            epochEEG = pop_loadset(epochEEGfilename);
            framesPerEpoch = epochEEG.pnts;
            totalTrials = epochEEG.trials;
            % load filtered EEG
            subStruct = load(filterEEGfilename);
            filteredStruct = subStruct.vgsBPfilterdata;
            thetaEEG = filteredStruct.vgsBPthetaEEG;
            lowbetaEEG = filteredStruct.vgsBPlowbetaEEG;
            highbetaEEG = filteredStruct.vgsBPhighbetaEEG;
            
            % average filtered time-courses across F-row
            thetaEEG_frow = squeeze(mean(thetaEEG([4 5 6 7 37 38 39 40 41],:),1));
            lowbetaEEG_frow = squeeze(mean(lowbetaEEG([4 5 6 7 37 38 39 40 41],:),1));
            highbetaEEG_frow = squeeze(mean(highbetaEEG([4 5 6 7 37 38 39 40 41],:),1));
            
            % reshape filtered EEG
            % original filtered EEG data is 1 x (totalTrials*framesPerEpoch)
            totalFrames = length(thetaEEG_frow);
            testNumTrials = totalFrames/framesPerEpoch;
            % check that computed number of trials matches actual number of trials
            if testNumTrials == totalTrials
                fprintf('Computed number of trials and total trials match\n')
            else
                fprintf('Computed number of trials and total trials do not match\n')
                noMatch(end+1,:) = [suID,scanDate];
                continue;
            end
            % reshape filtered EEG data to matrix of trials*times
            thetaEEG_frow_trials = reshape(thetaEEG_frow,[totalTrials,framesPerEpoch]);
            lowbetaEEG_frow_trials = reshape(lowbetaEEG_frow,[totalTrials,framesPerEpoch]);
            highbetaEEG_frow_trials = reshape(highbetaEEG_frow,[totalTrials,framesPerEpoch]);

            % truncate filtered EEGs to start (O ms) and end (500 ms) of prep period
            time0 = find(epochEEG.times==0);
            time500 = find(epochEEG.times==500);
            thetaEEG_frow_trials = thetaEEG_frow_trials(:,time0:time500);
            lowbetaEEG_frow_trials = lowbetaEEG_frow_trials(:,time0:time500);
            highbetaEEG_frow_trials = highbetaEEG_frow_trials(:,time0:time500);

            %% Spontaneous CFC- VGS Trials
            % sliding window and computing correlation in each window
            framesPerWindow = 15; % samples (100 ms)
            overlap = 10; % samples (66.67 ms)
            for currentTrial = 1:totalTrials
               % Generating time windows for each trial
                windowed_thetaEEG = buffer(thetaEEG_frow_trials(currentTrial,:),framesPerWindow,overlap,'nodelay');
                windowed_lowbetaEEG = buffer(lowbetaEEG_frow_trials(currentTrial,:),framesPerWindow,overlap,'nodelay');
                windowed_highbetaEEG = buffer(highbetaEEG_frow_trials(currentTrial,:),framesPerWindow,overlap,'nodelay');
                % buffer() outputs frames per window x # of windows
                % columns are the time windows, rows are the frames in the window
                numWindows = size(windowed_thetaEEG,2);
                
                % In each window, calculating spontaneous CFC for each trial
                for currentWindow = 1:numWindows
                    currentWindow_thetaEEG = nonzeros(windowed_thetaEEG(:,currentWindow));
                    currentWindow_lowbetaEEG = nonzeros(windowed_lowbetaEEG(:,currentWindow));
                    currentWindow_highbetaEEG = nonzeros(windowed_highbetaEEG(:,currentWindow));
                    mat = [currentWindow_thetaEEG,currentWindow_lowbetaEEG,currentWindow_highbetaEEG];
                    subvgsR_spontaneous(currentTrial,currentWindow,:,:) = corrcoef(mat);
                end
            end
            % Average subvgsR_spontaneous across trials to get average spontaneous correlation for each subject
            vgsR_spontaneous(currentEEG,:,:,:) = squeeze(mean(subvgsR_spontaneous,1));

            %% Evoked CFC- VGS Trials
            % average filtered EEGs across trials to get evoked correlation
            thetaEEG_frow_avgtrials = mean(thetaEEG_frow_trials,1);
            lowbetaEEG_frow_avgtrials = mean(lowbetaEEG_frow_trials,1);
            highbetaEEG_frow_avgtrials = mean(highbetaEEG_frow_trials,1);
            windowed_thetaEEG_avgtrials = buffer(thetaEEG_frow_avgtrials,framesPerWindow,overlap,'nodelay');
            windowed_lowbetaEEG_avgtrials = buffer(lowbetaEEG_frow_avgtrials,framesPerWindow,overlap,'nodelay');
            windowed_highbetaEEG_avgtrials = buffer(highbetaEEG_frow_avgtrials,framesPerWindow,overlap,'nodelay');
            % buffer() outputs frames per window x # of windows
            % columns are the time windows, rows are the frames in the window
            numWindows = size(windowed_thetaEEG_avgtrials,2);
            
            % evoked CFC at each time window
            for currentWindow = 1:numWindows
                % nonzeros() b/c buffer zero-pads windows if not evenly divisible
                currentWindow_thetaEEG_avgtrials = nonzeros(windowed_thetaEEG_avgtrials(:,currentWindow));
                currentWindow_lowbetaEEG_avgtrials = nonzeros(windowed_lowbetaEEG_avgtrials(:,currentWindow));
                currentWindow_highbetaEEG_avgtrials = nonzeros(windowed_highbetaEEG_avgtrials(:,currentWindow));
                mat = [currentWindow_thetaEEG_avgtrials,currentWindow_lowbetaEEG_avgtrials,currentWindow_highbetaEEG_avgtrials];
                vgsR_evoked(currentEEG,currentWindow,:,:) = corrcoef(mat);
            end
        else
            fprintf('Epoch EEG and filtered EEG do not match subject IDs; skipping')
        end
    end
    % save evoked correlation matrix
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgsCorrelationMats.mat','vgsR_evoked','vgsR_spontaneous')
end

%% Visualizing correlations
close all
times = linspace(50,479,14);
% Theta- Low Beta
corR_evoked_thetalowbeta = squeeze(corR_evoked(:,:,1,2));
vgsR_evoked_thetalowbeta = squeeze(vgsR_evoked(:,:,1,2));
corR_spontaneous_thetalowbeta = squeeze(corR_spontaneous(:,:,1,2));
vgsR_spontaneous_thetalowbeta = squeeze(vgsR_spontaneous(:,:,1,2));
figure
subplot(1,2,1)
errorbar(times,mean(corR_evoked_thetalowbeta,1),err_corevoked_thetalowbeta);
plot(times,corR_evoked_thetalowbeta,'LineWidth',1.5)
hold on
plot(times,vgsR_evoked_thetalowbeta,'LineWidth',1.5)
yline(0,'--k')
xlabel('Time (ms)'); ylabel('Correlation')
title('Evoked Coupling')
legend('Correct AS','VGS','Location','northwest')
ylim([-0.04 0.07]); grid on

subplot(1,2,2)
plot(times,corR_spontaneous_thetalowbeta,'LineWidth',1.5)
hold on
plot(times,vgsR_spontaneous_thetalowbeta,'LineWidth',1.5)
yline(0,'--k')
xlabel('Time (ms)'); ylabel('Correlation')
title('Spontaneous Coupling')
legend('Correct AS','VGS','Location','northwest')
ylim([-0.04 0.07]); grid on

sgtitle('Theta-Low Beta Coupling')

% Theta- High Beta
corR_evoked_thetahighbeta = mean(squeeze(corR_evoked(:,:,1,3)),1);
vgsR_evoked_thetahighbeta = mean(squeeze(vgsR_evoked(:,:,1,3)),1);
corR_spontaneous_thetahighbeta = mean(squeeze(corR_spontaneous(:,:,1,3)),1);
vgsR_spontaneous_thetahighbeta = mean(squeeze(vgsR_spontaneous(:,:,1,3)),1);
figure
subplot(1,2,1)
plot(times,corR_evoked_thetahighbeta,'LineWidth',1.5)
hold on
plot(times,vgsR_evoked_thetahighbeta,'LineWidth',1.5)
yline(0,'--k')
xlabel('Time (ms)'); ylabel('Correlation')
title('Evoked Coupling')
legend('Correct AS','VGS','Location','northwest')
ylim([-0.06 0.015]); grid on

subplot(1,2,2)
plot(times,corR_spontaneous_thetahighbeta,'LineWidth',1.5)
hold on
plot(times,vgsR_spontaneous_thetahighbeta,'LineWidth',1.5)
yline(0,'--k')
xlabel('Time (ms)'); ylabel('Correlation')
title('Spontaneous Coupling')
legend('Correct AS','VGS','Location','northwest')
ylim([-0.06 0.015]); grid on

sgtitle('Theta-High Beta Coupling')

% Low Beta- High Beta
corR_evoked_lowbetahighbeta = mean(squeeze(corR_evoked(:,:,2,3)),1);
vgsR_evoked_lowbetahighbeta = mean(squeeze(vgsR_evoked(:,:,2,3)),1);
corR_spontaneous_lowbetahighbeta = mean(squeeze(corR_spontaneous(:,:,2,3)),1);
vgsR_spontaneous_lowbetahighbeta = mean(squeeze(vgsR_spontaneous(:,:,2,3)),1);
figure
subplot(1,2,1)
plot(times,corR_evoked_lowbetahighbeta,'LineWidth',1.5)
hold on
plot(times,vgsR_evoked_lowbetahighbeta,'LineWidth',1.5)
yline(0,'--k')
xlabel('Time (ms)'); ylabel('Correlation')
title('Evoked Coupling')
legend('Correct AS','VGS','Location','southwest')
ylim([0 0.34]); grid on

subplot(1,2,2)
plot(times,corR_spontaneous_lowbetahighbeta,'LineWidth',1.5)
hold on
plot(times,vgsR_spontaneous_lowbetahighbeta,'LineWidth',1.5)
yline(0,'--k')
xlabel('Time (ms)'); ylabel('Correlation')
title('Spontaneous Coupling')
legend('Correct AS','VGS','Location','southwest')
ylim([0 0.34]); grid on

sgtitle('Low Beta-High Beta Coupling')