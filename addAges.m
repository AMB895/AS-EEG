%% Adding age and sex to preprocessed EEG structures
clear
% adding paths
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable_20250218.mat') % to get ages
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
eeglab;

% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
% datafolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
datafolder = [taskdirectory,'/AfterWhole/ICAwholeClean_homogenize/'];
% dataname = '*.set';
dataname = '*icapru.set';

EEGfilenames = dir([datafolder,dataname]);
for currentEEG=1:size(EEGfilenames)
    
    % defining input EEG file from epochclean_homogenize
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [datafolder,filename];
    
    % getting subID and scan date
    [~, currentName, ~ ] = fileparts(inputfile);
    splitCurrentName = split(currentName,'_');
    subID = cell2mat(splitCurrentName(1));
    scanDate = cell2mat(splitCurrentName(2));

    % finding age from ErrorLatency Table
    index = find(ErrorLatencyTable.("LunaID")==str2double(subID) & ErrorLatencyTable.("ScanDate")==str2double(scanDate));
    if isempty(index)
        fprintf('No age for %s %s\n',subID,scanDate)
        continue;
    end
    AGE = ErrorLatencyTable.Age(index);
    SEX = ErrorLatencyTable.Sex(index);

    % Loading currentEEG
    EEG = pop_loadset(inputfile);
    EEG.age = AGE;
    EEG.sex = SEX;

    % Saving currentEEG with age in EEG structure
    EEG = pop_saveset(EEG,'savemode','resave');
    clear AGE SEX
end



