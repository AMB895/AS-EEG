%% Generate Scalp Maps of Prep Period
% adding path to eeglab and starting eeglab
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
eeglab;
%%
close all;
% Time-Frequency Points to plot scalp maps
timeFreqs = zeros(4,2,4);
timeFreqs(:,:,1) = [50 5;100 5;250 5; 400 5];
timeFreqs(:,:,2) = [50 10;100 10;250 10; 400 10];
timeFreqs(:,:,3) = [50 15;100 15;250 15; 400 15];
timeFreqs(:,:,4) = [50 25;100 25;250 25; 400 25];
timeFreqsLabels = ["Theta","Alpha","Low Beta","High Beta"];

%% Correct AS Trials
% loading in channel locations
EEGfiles = dir('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/ICAwholeClean_homogenize/*.set');
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/ICAwholeClean_homogenize/')
currentEEG = EEGfiles(1).name;
EEG = pop_loadset(currentEEG);
titstr = 'Correct Trials';
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/corERSPdata.mat')
% tftopo() needs ersp data as [time,freq,chans,subjects]
erspdata = permute(corerspdata,[3 4 2 1]);
for i=1:size(timeFreqs,3)
    figure;
    tftopo(erspdata,times,freqs,'chanlocs',EEG.chanlocs,'timefreqs',squeeze(timeFreqs(:,:,i)),'mode','ave');
    sgtitle(sprintf("%s: %s",titstr,timeFreqsLabels(i)))
end

%% Error AS Trials
% loading in channel locations
EEGfiles = dir('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/ICAwholeClean_homogenize/*.set');
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/ICAwholeClean_homogenize/')
currentEEG = EEGfiles(1).name;
EEG = pop_loadset(currentEEG);
titstr = 'Error Trials';
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/errcorERSPdata.mat')
% tftopo() needs ersp data as [time,freq,chans,subjects]
erspdata = permute(errcorerspdata,[3 4 2 1]);
for i=1:size(timeFreqs,3)
    figure;
    tftopo(erspdata,times,freqs,'chanlocs',EEG.chanlocs,'timefreqs',squeeze(timeFreqs(:,:,i)),'mode','ave');
    sgtitle(sprintf("%s: %s",titstr,timeFreqsLabels(i)))
end

%% VGS Trials
% loading in channel locations
EEGfiles = dir('/Volumes/Hera/Abby/preprocessed_data/vgs/AfterWhole/ICAwholeClean_homogenize/*.set');
addpath('/Volumes/Hera/Abby/preprocessed_data/vgs/AfterWhole/ICAwholeClean_homogenize/')
currentEEG = EEGfiles(1).name;
EEG = pop_loadset(currentEEG);
titstr = 'VGS';
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs/vgsERSPdata.mat')
% tftopo() needs ersp data as [time,freq,chans,subjects]
erspdata = permute(vgserspdata,[3 4 2 1]);
for i=1:size(timeFreqs,3)
    figure;
    tftopo(erspdata,times,freqs,'chanlocs',EEG.chanlocs,'timefreqs',squeeze(timeFreqs(:,:,i)),'mode','ave');
    sgtitle(sprintf("%s: %s",titstr,timeFreqsLabels(i)))
end