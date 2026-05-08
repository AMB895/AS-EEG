%% ERSP and ITC for correct AS, VGS, and error corrected AS
% amb895 05/07/2026
clear; 

% add paths
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/AS_EEG/')
addpath('/Volumes/Hera/Abby/Resources/slanCM/')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab("nogui");

%% Set up directories/ load in relevant data
% Getting age and sex from merge7t
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')

% Load merged 7t EEG and AS score information
merged7t_eeg = readtable('/Volumes/Hera/Abby/AS_EEG/Results/merged7t_ASscore_table.csv');

% Main directory of preprocessed data
maindir = hera('Abby/preprocessed_data');

% loading template of channel labels
load('/Volumes/Hera/Abby/AS_EEG/templatechannellabels.mat')

% RERUN; 0 = do not rerun, 1 = rerun
RERUN = 0;

% stim or prep onset
eventtype = 'prep';

% Directories to AS and VGS preprocessed data
taskdirAS = [maindir '/' 'anti']; 
taskdirVGS = [maindir '/' 'vgs']; 

% Directories to AS and VGS epoched data
epochdirAS = [taskdirAS,'/AfterWhole/kept_epochclean_prep/'];
epochdirVGS = [taskdirVGS,'/AfterWhole/kept_epochclean_prep/'];

% Epoch data file names
epochcorAS_filename = ['*_epochs_kept_cor_' eventtype '.set'];
epocherrcorAS_filename = ['*_epochs_kept_errcor_' eventtype '.set'];
epochVGS_filename = ['*_epochs_kept_' eventtype '.set'];

% Add all eeg files to one big structure to loop over
EEGfilenames = {dir([epochdirAS,epochcorAS_filename]) dir([epochdirVGS,epochVGS_filename]) dir([epochdirAS,epocherrcorAS_filename]) };
tasks = ["corAS","VGS","errcorAS"];

% Baseline start and end/ event start and end times (in ms)
basestart = -700;
baseend = -400;
eventstart = 0;
eventend = 500;

% comptue average trial baseline and single trial baseline
% base_norm = 'avgtrial_ers';
% basenormfolder = 'avgtrl';
base_norm = 'singletrial';
basenormfolder = 'singletrl';

%% Loop through each task and subject
for t = 1:length(tasks)
    % Defining current task EEG
    taskfilenames = EEGfilenames{t};
    tasktotalEEGs = length(taskfilenames);
    currtask = tasks(t);
    parfor currentEEG = 1:tasktotalEEGs
        % getting subject's id and scan date
        filename = [taskfilenames(currentEEG).name];
        filepath = [taskfilenames(currentEEG).folder];
        splitfilename = strsplit(filename,'_');
        lunaid = cell2mat(splitfilename(1));
        eeg_date = cell2mat(splitfilename(2));
        % define file name to save subject's data
        subtasksavepath = sprintf('/Volumes/Hera/Abby/AS_EEG/ERSPdata/%s_data/%s/SubjectData/%s',currtask,eventtype,basenormfolder);
        subtasksavename = sprintf('%s_%s_%s_%s_ersp.mat',lunaid,eeg_date,currtask,eventtype);
        % if already have ERSP for subject skip
        if exist(fullfile(subtasksavepath,subtasksavename),'file') && RERUN ==0
            fprintf('skipping; computed ERSP and time/freq for %s_%s\n',lunaid,eeg_date)
        else
            % defining subject's merged7t table
            subIdx = find(merged7t_eeg.lunaid==str2double(lunaid) & merged7t_eeg.eeg_date==str2double(eeg_date));
            subTable = merged7t_eeg(subIdx,:);
            
            % skip if subject is not in merged 7t
            if isempty(subTable)
                warning('Subject %s %s is not in merge 7t; skipping\n',lunaid,eeg_date)
                continue
            end
            % skip if age is empty in merged 7t table
            if isnan(subTable.eeg_age)
                warning('Subject %s %s does not have age; skipping\n',lunaid,eeg_date)
                continue
            end
            % load subject's EEG
            EEG = pop_loadset('filename',filename,'filepath',filepath);
            chanNames = string({EEG.chanlocs.labels}');
            % skipping if channel labels do not match template
            if all(tempChanLabels~=chanNames)
                warning('Check Subject Channel Labels %s %s; skipping\n',lunaid,eeg_date)
                continue
            end
            % skip if subject has 0 or 1 epochs (fails if subject has only 1 epoch)
            numEpochs = length(EEG.epoch);
            if numEpochs < 1
                warning('Subject %s %s has %d epochs; skipping\n',lunaid,eeg_date,numEpochs)
                continue
            end
            
            % define when event starts and ends to reduce size of 4d array for saving
            eventstart_idx = find(EEG.times>=eventstart,1,'first');
            eventend_idx = find(EEG.times>=eventend,1,'first');
            
            % Compute TF reconstruction of power
            [sub_totpower, times, freqs, sub_itc] = F_tfpower(EEG,'itc',1);
            % Comptue ERSP from power
            [sub_ersp, ersptimes] = F_ersp(sub_totpower,times,[basestart baseend],[eventstart eventend],'base_norm',base_norm);
            
            % Just prep period of power and itc
            sub_totpower = sub_totpower(:,eventstart_idx:eventend_idx,:,:);
            sub_itc = sub_itc(:,eventstart_idx:eventend_idx,:);
            
            % Save subject data
            sub_struct = struct('sub_ersp',sub_ersp,...
                'sub_itc',sub_itc, 'sub_totpower', sub_totpower,...
                'times',ersptimes,'freqs',freqs);
            parsave(fullfile(subtasksavepath,subtasksavename),sub_struct)            
        end
        % clear variables for next subject
        sub_totpower = []; sub_ersp = []; sub_itc = [];
        ersptimes = []; times = []; freqs = [];
        
    end
    
end

%% function to save within parfor loop
function parsave(filename, sub_struct)
save(filename,'-fromstruct',sub_struct);
end