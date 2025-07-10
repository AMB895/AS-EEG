%% Adding age and sex to preprocessed EEG structures
clear
% adding paths
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
eeglab;

% Main directory
maindir = hera('Abby/preprocessed_data');
% task = 'anti';
task = 'vgs';

if task == "anti"
    ageTable = load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable_20250320.mat'); % to get ages
    ageTable = ageTable.ErrorLatencyTable;
elseif task == "vgs"
    addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/') % for merge 7t
    merged7t=readtable('merged_7t.csv'); % to get age and sex
    
    % only care about lunaid, eeg scan date, age at eeg scan and sex from merged7t
    merged7t_eeg = [merged7t(:,'lunaid'), merged7t(:,'eeg_date'), merged7t(:,'eeg_age'), merged7t(:,'sex')];
    
    % finding subjects who did not do not have an eeg scan and removing them from merge_7t_eeg
    noeeg_idx = find(isnan(merged7t_eeg.eeg_date));
    merged7t_eeg(noeeg_idx,:) = [];

    merged7t_eeg = renamevars(merged7t_eeg,["lunaid","eeg_date","eeg_age","sex"],["LunaID","ScanDate","Age","Sex"]);
    % changing column names 
    ageTable = merged7t_eeg;
end

taskdirectory = [maindir '/' task]; 
% In path
datafolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
% datafolder = [taskdirectory,'/AfterWhole/ICAwholeClean_homogenize/'];
dataname = '*.set';
% dataname = '*icapru.set';

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
    index = find(ageTable.("LunaID")==str2double(subID) & ageTable.("ScanDate")==str2double(scanDate));
    if isempty(index)
        fprintf('No age for %s %s\n',subID,scanDate)
        continue;
    end
    AGE = ageTable.Age(index);
    SEX = ageTable.Sex(index);

    % Loading currentEEG
    EEG = pop_loadset(inputfile);
    EEG.age = AGE;
    EEG.sex = SEX;

    % Saving currentEEG with age in EEG structure
    EEG = pop_saveset(EEG,'savemode','resave');
    clear AGE SEX
end



