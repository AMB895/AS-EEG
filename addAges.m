%% Adding ages to preprocessed EEG structures for Study analysis
% adding paths
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat') % to get ages
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')


% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
epochedFolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedName = '*_epochs_kept.set';

EEGfilenames = dir([epochedFolder,epochedName]);
for currentEEG=1:size(EEGfilenames)
    
    % defining input EEG file from epochclean_homogenize
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedFolder,filename];
    
    % getting subID and scan date
    [~, currentName, ~ ] = fileparts(inputfile);
    splitCurrentName = split(currentName,'_');
    subID = cell2mat(splitCurrentName(1));
    scanDate = cell2mat(splitCurrentName(2));

    % finding age from ErrorLatency Table
    index = find(ErrorLatencyTable.("Luna ID")==str2double(subID) & ErrorLatencyTable.("Scan Date")==str2double(scanDate));
    if isempty(index)
        fprintf('No age for %s %s\n',subID,scanDate)
        continue;
    end
    AGE = ErrorLatencyTable.Age(index);

    % Loading currentEEG
    EEG = pop_loadset(inputfile);
    EEG.age = AGE;

    % Saving currentEEG with age in EEG structure
    EEG = pop_saveset(EEG,'savemode','resave');
    clear AGE
end


% Creating a new Study with EEGLAB


