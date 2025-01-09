%% Threshold-Free Cluster Enhancement Method
% From limo_tools
% Functions:
% channeighbstructmat = limo_neighbourdist(EEG,neighbourdist)
% ^ creates neighbourhood distance matrix to control for mulitple comparisions
%   neighbourdist: a number that indicates the neighbour distance for biosemi 64
%       ^ fieldtrip may have the neighbour matrix already
% [tfce_score, thresholded_maps] = limo_tfce(type, data,channeighbstructmat, updatebar, E, H, dh)
%   type: 1 for 1D data, 2 for 2D data, 3 for 3D data
%   data: map of t/F values or set of t/F maps computed under H0
%   channeighbstructmat: neighbourhood matrix for clustering (size is # of chans x # of chans)

% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/STUDY/') % to get ERSP,ITC & ERP data
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

% loading data
load('erspSTUDY.mat')
load('itcSTUDY.mat')
ftNeighbours = load('biosemi64_neighb.mat');

% STUDY file path and name
studypath = '/Volumes/Hera/Abby/AS_EEG/STUDY/';
studyname = 'AS_correct_PrepPeriod.study';

% % load study
% [STUDY, ALLEEG] = pop_loadstudy('filename',studyname,'filepath',studypath);

% get EEGfilenames
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
epochedFolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedName = '*_epochs_kept.set';
EEGfilenames = dir([epochedFolder,epochedName]);

%% Make neighbourhood distance matrix for each subject

for currentEEG = 1:size(EEGfilenames)
    % defining input EEG file from epochclean_homogenize
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedFolder,filename];
    
    % getting subID and scan date
    [~, currentName, ~ ] = fileparts(inputfile);
    splitCurrentName = split(currentName,'_');
    subID = cell2mat(splitCurrentName(1));
    scanDate = cell2mat(splitCurrentName(2));

    % loading subject's EEG
    EEG = pop_loadset(inputfile);

    % trying to get channel neighbours
    eeglabpath = fileparts(which('eeglab'));
    cap_location = fullfile(eeglabpath,'/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
    EEG=pop_chanedit(EEG, 'lookup', cap_location);
end

%% trying to get channel neighbours
EEG = pop_loadset('10129_20180919_anti_Rem.set','/Volumes/Hera/Abby/preprocessed_data/anti/remarked/');
% trying to get channel neighbours
    eeglabpath = fileparts(which('eeglab'));
    cap_location = fullfile(eeglabpath,'/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
    EEG=pop_chanedit(EEG, 'lookup', cap_location);
    pop_chanedit(EEG)
channelLabels = {EEG.chanlocs(1:64).labels};
save('/Volumes/Hera/Abby/AS_EEG/STUDY/channelLabels.mat',"channelLabels")