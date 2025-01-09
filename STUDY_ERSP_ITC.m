%% ERSP/ ITC Study
% add paths
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

% epoched data path
epochpath = '/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/';

% STUDY file path and name
studypath = '/Volumes/Hera/Abby/AS_EEG/STUDY/';
studyname = 'AS_correct_PrepPeriod.study';

% load study
[STUDY, ALLEEG] = pop_loadstudy('filename',studyname,'filepath',studypath);

%% make sure all channel labels are the same across subjects
load('channelLabels.mat');
% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
epochedFolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedName = '*_epochs_kept.set';
EEGfilenames = dir([epochedFolder,epochedName]);

for currentEEG = 1:size(EEGfilenames)
    % defining current EEG
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedFolder,filename];
    [d, currentName, ext ] = fileparts(inputfile);
    % loading subject's EEG
    EEG = pop_loadset(inputfile);
    numChans = EEG.nbchan;
    for chan = 1:numChans
        EEG.chanlocs(chan).labels = channelLabels{chan};
        EEG.chanlocs(chan).numLabels = chan;
    end
    EEG = pop_saveset(EEG,'filename',filename,'filepath',epochedFolder);
end


%% Precompute ERSP & ITC, ERP, and Spec
% ERSP & ITC
[STUDY,ALLEEG] = std_precomp(STUDY,ALLEEG,'channels','ersp','on','itc','on','erspparams',{...
    'cycles',[3 0.5],'freqs',[3 50],'baseline',[-700 -400]});
% save study
STUDY = pop_savestudy(STUDY,[],'filename',studyname,'filepath',studypath);

% ERP & Spec
[STUDY,ALLEEG] = std_precomp(STUDY,ALLEEG,'channels','erp','on','erpparams',{'rmbase',[-700 -400]},...
    'spec','on');
% save study
STUDY = pop_savestudy(STUDY,[],'filename',studyname,'filepath',studypath);
%% delete files
% erspdata = dir([epochpath,'*.dattimef']);
% erpdata = dir([epochpath,'*.daterp']);
% specdata = dir([epochpath,'*.datspec']);
% cd /Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/
% delete(erspdata.name)
% delete(erpdata.name)
% delete(specdata.name)
%% Read ERSP and ITC data from STUDY
cd /Volumes/Hera/Abby/AS_EEG/STUDY/
% checking to see if data has already been saved
if exist('erspSTUDY.mat','file')
    fprintf('ERSP data already saved\n')
else
    % erspdata 4D matrix [Freqs,Times,Channels,Subjects]
    [STUDY,erspdata,timesERSP,freqsERSP] = std_readdata(STUDY,ALLEEG,'channels',{ALLEEG(1).chanlocs.labels},'datatype','ersp');
    save(fullfile(studypath,'erspSTUDY.mat'),'erspdata','timesERSP','freqsERSP')
end

if exist('itcSTUDY.mat','file')
    fprintf('ITC data already saved\n')
else
    % itcdata 4D matrix [Freq, Time, Channel, Subject]
    [STUDY,itcdata,timesITC,freqsITC] = std_readdata(STUDY,ALLEEG,'channels',{ALLEEG(1).chanlocs.labels},'datatype','itc');
    save(fullfile(studypath,'itcSTUDY.mat'),'itcdata','timesITC','freqsITC')
end

if exist('erpSTUDY.mat','file')
    fprintf('ERP data already saved\n')
else
    % erpdata 3D matrix [Time,Channel,Subject]
    [STUDY,erpdata,timesERP] = std_readdata(STUDY,ALLEEG,'channels',{ALLEEG(1).chanlocs.labels},'datatype','erp');
    save(fullfile(studypath,'erpSTUDY.mat'),'erpdata','timesERP')
end

if exist('specSTUDY.mat','file')
    fprintf('Spec data already saved\n')
else
    % specdata 3D matrix [Freq, Channel, Subject]
    % [STUDY,specdata,freqsSPEC] = std_readdata(STUDY,ALLEEG,'channels',{ALLEEG(1).chanlocs.labels},'datatype','spec');
    % Error with 11651_20190816.datspec File might be corrupt
    % save(fullfile(studypath,'specSTUDY.mat'),'specdata','freqsSPEC')
end

