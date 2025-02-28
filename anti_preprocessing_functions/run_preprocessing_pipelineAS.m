function [] = run_preprocessing_pipelineAS(task)

% load in paths to feildtrip and eeglab
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/'));
addpath(genpath('/resources/Euge/'))
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath(genpath('/Volumes/Hera/Abby/AS_EEG/anti_preprocessing_functions/'))
% addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20220104/')
%run('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/resources/eeglab2022.1/eeglab.m');
run('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/eeglab.m')

% Outpath
%maindir = hera('Projects/7TBrainMech/scripts/eeg/Shane/preprocessed_data');
maindir = hera('Abby/preprocessed_data');
task = 'anti';

disp(["i am running" task])

taskdirectory = [maindir '/' task];

% initial values
lowBP = 0.5;
topBP = 70;
FLAG = 1;

%% settings
only128 = 0; % 0==do all, 1==only 128 channel subjects
condition = 1; %0 - if you want to overwrite an already existing file; 1- if you want it to skip subjects who have already been run
dryrun = 0;

% remarking antisaccade trials
% chages trials to be a single digit
remarkAS(taskdirectory,dryrun)

% gather all file paths for all subjects remarked data
setfilesDir = [taskdirectory, '/remarked/*.set'];
setfiles = all_remarked_setAS(setfilesDir);
n = length(setfiles);%number of EEG sets to preprocessa

% saving input files that failed to be preprocessed
errorProcessing=[];

% loop through every subject to preprocess
for i = 1:n
    inputfile = setfiles{i};
    try
        preprocessing_pipelineAS(inputfile,taskdirectory,lowBP,topBP,FLAG,condition,task)
    catch e
        fprintf('Error processing "%s": %s\n',inputfile, e.message)
        errorProcessing(end+1,:) = [inputfile];
        for s=e.stack
            disp(s)
        end
    end
end


%% select ICA values to reject

ICA_Path = [taskdirectory '/ICAwhole'];
CleanICApath = [taskdirectory '/AfterWhole/ICAwholeClean/'];


EEGfileNames = dir([ICA_Path '/*.set']);

for fidx = 1:length(EEGfileNames)
    filename = EEGfileNames(fidx).name;
    locs = file_locs(fullfile(ICA_Path,filename), taskdirectory, task);
    if exist(locs.ICAwholeClean, 'file')
        fprintf('skipping; already created %s\n', locs.ICAwholeClean);
        continue
    end

    selectcompICA
end


%% Homogenize Chanloc
datapath = [taskdirectory '/AfterWhole/ICAwholeClean'];
savepath = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize'];

setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};
for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

correction_cap_location = hera('Projects/7TBrainMech/scripts/eeg/Shane/resources/ELchanLoc.ced');
for i = 1:length(setfiles)
    homogenizeChanLoc(setfiles{i},correction_cap_location,savepath, taskdirectory, task)
end


%% filter out the 60hz artifact from electronics
datapath = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize'];

setfiles0 = dir([datapath,'/*icapru.set']);
setfiles = {};
for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
end

for i = 1:length(setfiles)
    inputfile = setfiles{i};
    [filepath,filename ,ext] =  fileparts((setfiles{i}));

    EEG = pop_loadset(inputfile);
    EEG = pop_eegfiltnew(EEG, 59, 61, [], 1, [], 0);
    EEG = pop_saveset(EEG, 'filename', filename, 'filepath', datapath);
end

%% Clean epochs and remove ones that are bad
path_data = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize/'];
epoch_path = [taskdirectory '/AfterWhole/ICAwholeClean_homogenize/'];
epoch_folder = [taskdirectory '/AfterWhole/epoch/'];
kept_epoch_folder = [taskdirectory '/AfterWhole/epochclean_homogenize/'];

EEGfileNames = dir([path_data, '/*_icapru.set']);

revisar = {};
% Define what type of trial you want before running!!!
% trialtype = 1; % correct trials
% trialtype = 0; % incorrect trials
trialtype = 2; % error corrected trials

for currentEEG=1:size(EEGfileNames,1)
    filename=[EEGfileNames(currentEEG).name];
    inputfile = [epoch_path,filename];
    try
        revisar{currentEEG} = epochcleanAS(inputfile,epoch_folder,kept_epoch_folder,task,taskdirectory,trialtype);
    catch e
        fprintf('Error processing "%s": %s\n',inputfile, e.message)
        for s=e.stack
            disp(s)
        end
    end
end