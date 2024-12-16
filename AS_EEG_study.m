%% Creating a new Study with EEGLAB
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

% Out path
studyFolder = '/Volumes/Hera/Abby/AS_EEG/STUDY/';

% Files being created
alleegName = 'ALLEEG_study.mat';

% First need to group all subject's EEG structures into one ALLEEG and save in STUDY folder
EEGfilenames = dir([epochedFolder,epochedName]);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
for currentEEG = 1:size(EEGfilenames)
    % defining input EEG file from epochclean_homogenize
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedFolder,filename];
    
    % getting subID and scan date
    [~, currentName, ~ ] = fileparts(inputfile);
    splitCurrentName = split(currentName,'_');
    subID = cell2mat(splitCurrentName(1));
    scanDate = cell2mat(splitCurrentName(2));
    subject = sprintf('%s_%s',subID,scanDate);

    % loading currentEEG
    EEG = pop_loadset(inputfile);

    % Don't add files with 0 epochs
    numOfepochs = length(EEG.epoch);
    if numOfepochs<1
        fprintf('skipping; No epochs %s %s',subID,scanDate)
        continue;
    end

    % Adding currentEEG to ALLEEG
    [ALLEEG ,EEG ,CURRENTSET] = eeg_store(ALLEEG,EEG);

    % Adding in empty variables to ALLEEG
    if isempty(ALLEEG(CURRENTSET).filename)
        ALLEEG(CURRENTSET).filename = filename;
    end

    if isempty(ALLEEG(CURRENTSET).filepath)
        ALLEEG(CURRENTSET).filepath = inputfile;
    end

    if isempty(ALLEEG(CURRENTSET).subject)
        ALLEEG(CURRENTSET).subject = subject;
    end

    if isempty(ALLEEG(CURRENTSET).condition)
        ALLEEG(CURRENTSET).condition = 'anti_Rem';
    end

end
save(fullfile(studyFolder,alleegName),'ALLEEG')
%% Making STUDY
[STUDY] = std_makedesign(STUDY,ALLEEG);