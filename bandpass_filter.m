%% Filtering epoched EEG data in significant frequency clusters from ERSP/TFCE
clear;close all;
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/preprocessed_data/vgs/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ClusterStats/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;
% load in mask of significant activation in correct AS trials
load('corGroupActClusters.mat')
load('CorrectTrials/corERSPdata.mat','freqs')
%% getting upper and lower frequencies for each cluster
% number the clusters
groupact_clusters_cor = mask_cor_groupact.*t_cor_groupact;
[labeled_clusters, ~] = bwlabel(groupact_clusters_cor);

% bwlabel() combines low beta and high beta cluster
% low beta cluster high freq cutoff = 22.53 Hz
freq22 = find(freqs>22 & freqs<23,1);
[idxHighBetafreq,idxHighBetatime] = find(labeled_clusters==2);
% making the high beta cluster separate from low beta cluster
for currentTime = 1:length(idxHighBetatime)
for currentFreq = 1:length(idxHighBetafreq)
   if idxHighBetafreq(currentFreq)>=freq22 && labeled_clusters(idxHighBetafreq(currentFreq),idxHighBetatime(currentTime))==2
       labeled_clusters(idxHighBetafreq(currentFreq),idxHighBetatime(currentTime)) = 3;
   end
end
end
numClusters = 3;
% finding upper and lower frequencies for each cluster
for currentCluster = 1:numClusters
    % find individual clusters
    cluster_mask_cor(currentCluster,:,:) = labeled_clusters == currentCluster;
    % determine frequency band
    [freqIdx,timeIdx] = find(squeeze(cluster_mask_cor(currentCluster,:,:)));
    minFreqIdx = min(freqIdx);
    maxFreqIdx = max(freqIdx);
    minFreq = freqs(minFreqIdx);
    maxFreq = freqs(maxFreqIdx);
    % save frequency bands and indicies
    clusterFreqbands(currentCluster,1) = minFreq;
    clusterFreqbands(currentCluster,2) = maxFreq;
end

%% Load in each subject's EEG and bandpass filter Correct AS Trials
% Main directory of epoched EEG data
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% Epoched data path and name
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept_cor.set';
% Save filtered EEG path
filteredEEGpath = '/Volumes/Hera/Abby/AS_EEG/BPfilteredData/';
% Getting all epoched data files
EEGfilenames = dir([epochedpath,epochedname]);
totalEEGs = size(EEGfilenames,1);

for currentEEG = 1:length(EEGfilenames)
    % getting file name of EEG and subect's ID and scan date
    eegfilename = [EEGfilenames(currentEEG).name];
    splitfilename = strsplit(eegfilename,'_');
    subID = cell2mat(splitfilename(1));
    scanDate = cell2mat(splitfilename(2));
    % loading subject's EEG
    EEG = pop_loadset(eegfilename);
    % save names
    filterdatasavename = [subID,'_',scanDate,'_','corBPfilteredEEG.mat'];
    
    if exist(fullfile(filteredEEGpath,filterdatasavename),'file')
        fprintf('skipping; Already band-pass filtered subject %s %s EEG\n',subID,scanDate)
    else
        % BP filtering theta band
        thetaLF = clusterFreqbands(1,1); % eegfilt() fails when low freq is 3 but works when low freq is 3.001????
        thetaHF = clusterFreqbands(1,2);
        % filter data from all channels
        corthetaEEG = eegfilt(EEG.data(:,:,:),EEG.srate,3.001,thetaHF,EEG.pnts,[],0,'fir1',0);
    
        % BP filterig low beta band
        lowbetaLF = clusterFreqbands(2,1);
        lowbetaHF = clusterFreqbands(2,2);
        % filter data from all channels
        corlowbetaEEG = eegfilt(EEG.data(:,:,:),EEG.srate,lowbetaLF,lowbetaHF,EEG.pnts,[],0,'fir1',0);
    
        % BP filtering high beta band
        highbetaLF = clusterFreqbands(3,1);
        highbetaHF = clusterFreqbands(3,2);
        % filter data from all channels
        corhighbetaEEG = eegfilt(EEG.data(:,:,:), EEG.srate,highbetaLF,highbetaHF,EEG.pnts,[],0,'fir1',0);
        
        corBPfilterdata.corBPthetaEEG = corthetaEEG;
        corBPfilterdata.corBPlowbetaEEG = corlowbetaEEG;
        corBPfilterdata.corBPhighbetaEEG = corhighbetaEEG;
        save(fullfile(filteredEEGpath,filterdatasavename),'corBPfilterdata')
        
        clear corBPfilterdata
    end
end

%% Load in each subject's EEG and bandpass filter VGS Trials
% Main directory of epoched EEG data
maindir = hera('Abby/preprocessed_data');
task = 'vgs';
taskdirectory = [maindir '/' task]; 
% Epoched data path and name
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept.set';
% Save filtered EEG path
filteredEEGpath = '/Volumes/Hera/Abby/AS_EEG/BPfilteredData/';
% Getting all epoched data files
EEGfilenames = dir([epochedpath,epochedname]);
totalEEGs = size(EEGfilenames,1);

for currentEEG = 1:length(EEGfilenames)
    % getting file name of EEG and subect's ID and scan date
    eegfilename = [EEGfilenames(currentEEG).name];
    splitfilename = strsplit(eegfilename,'_');
    subID = cell2mat(splitfilename(1));
    scanDate = cell2mat(splitfilename(2));
    % loading subject's EEG
    EEG = pop_loadset(eegfilename);
    % save names
    filterdatasavename = [subID,'_',scanDate,'_','vgsBPfilteredEEG.mat'];
    
    if exist(fullfile(filteredEEGpath,filterdatasavename),'file')
        fprintf('skipping; Already band-pass filtered subject %s %s EEG\n',subID,scanDate)
    else
        % BP filtering theta band
        thetaLF = clusterFreqbands(1,1); % eegfilt() fails when low freq is 3 but works when low freq is 3.001????
        thetaHF = clusterFreqbands(1,2);
        % filter data from all channels
        vgsthetaEEG = eegfilt(EEG.data(:,:,:),EEG.srate,3.001,thetaHF,EEG.pnts,[],0,'fir1',0);
    
        % BP filterig low beta band
        lowbetaLF = clusterFreqbands(2,1);
        lowbetaHF = clusterFreqbands(2,2);
        % filter data from all channels
        vgslowbetaEEG = eegfilt(EEG.data(:,:,:),EEG.srate,lowbetaLF,lowbetaHF,EEG.pnts,[],0,'fir1',0);
    
        % BP filtering high beta band
        highbetaLF = clusterFreqbands(3,1);
        highbetaHF = clusterFreqbands(3,2);
        % filter data from all channels
        vgshighbetaEEG = eegfilt(EEG.data(:,:,:), EEG.srate,highbetaLF,highbetaHF,EEG.pnts,[],0,'fir1',0);
        
        vgsBPfilterdata.vgsBPthetaEEG = vgsthetaEEG;
        vgsBPfilterdata.vgsBPlowbetaEEG = vgslowbetaEEG;
        vgsBPfilterdata.vgsBPhighbetaEEG = vgshighbetaEEG;
        save(fullfile(filteredEEGpath,filterdatasavename),'vgsBPfilterdata')
        
        clear vgsBPfilterdata
    end
end