%% Generate and compile spontaneous and evoked power
% amb895 05/08/2026
clear;

% adding necessary paths 
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/AS_EEG/')

%% Setting up directories and settings
% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab("nogui");

% Getting age and sex from merge7t
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
resultsdir = hera('Abby/AS_EEG/Results');
% load merged7t AS score table
load(fullfile(resultsdir,'merged7t_ASscore_table.mat'))  

% Main directory of preprocessed data
maindirpreproc = hera('Abby/preprocessed_data');

% loading template of channel labels
load('/Volumes/Hera/Abby/AS_EEG/templatechannellabels.mat')

% RERUN; 0 = do not rerun, 1 = rerun
RERUN = 1;

frow = [4 5 6 7 37 38 39 40 41];

% stim or prep onset
eventtype = 'prep';

% directories to AS and VGS preprocessed data
taskdirAS = [maindirpreproc '/' 'anti']; 
taskdirVGS = [maindirpreproc '/' 'vgs'];

% directories to AS and VGS epoched data
epochdirAS = [taskdirAS '/AfterWhole/kept_epochclean_prep/'];
epochdirVGS = [taskdirVGS '/AfterWhole/kept_epochclean_prep/'];

% epoched data name for correct AS and VGS
epochcorAS_filename = ['*_epochs_kept_cor_' eventtype '.set'];
epochVGS_filename = ['*_epochs_kept_' eventtype '.set'];

% large cell of all file names
EEGfilenames = {dir([epochdirAS,epochcorAS_filename]), dir([epochdirVGS, epochVGS_filename])};

% baseline start/ end and event start/end
basestart = -700;
baseend = -400;
eventstart = 0;
eventend = 500;

% define tasks
tasks = ["corAS", "VGS"];

%% loop over tasks and subjects
temptable = merged7t_ASscore_table;
N = height(temptable);
sponpower = NaN(N,1);
evopower = NaN(N,1);
totpower = NaN(N,1);
for t = 1:length(tasks)
    taskfilenames = EEGfilenames{t};
    tasktotalEEGs = length(EEGfilenames{t}); 
    currtask = char(tasks(t));
    for currentEEG = 1:tasktotalEEGs
        % get subject's id and eeg date
        filename = [taskfilenames(currentEEG).name];
        filepath = [taskfilenames(currentEEG).folder];
        splitfilename = strsplit(filename,'_');
        lunaid = cell2mat(splitfilename(1));
        eeg_date = cell2mat(splitfilename(2));
        % defining subject's merged7t table
        subIdx = find(merged7t_ASscore_table.lunaid==str2double(lunaid) & merged7t_ASscore_table.eeg_date==str2double(eeg_date));
        subTable = merged7t_ASscore_table(subIdx,:);

        % skip if subject is not in merged 7t
        if isempty(subTable)
            warning('Subject %s %s is not in merge 7t; skipping\n',lunaid,eeg_date)
            sponpower(subIdx,1) = NaN;
            evopower(subIdx,1) = NaN;
            totpower(subIdx,1) = NaN;
            continue
        end
        % skip if age is empty in merged 7t table
        if isnan(subTable.eeg_age)
            warning('Subject %s %s does not have age; skipping\n',lunaid,eeg_date)
            sponpower(subIdx,1) = NaN;
            evopower(subIdx,1) = NaN;
            totpower(subIdx,1) = NaN;
            continue
        end
        % load subject's EEG
        EEG = pop_loadset('filename',filename,'filepath',filepath);
        chanNames = string({EEG.chanlocs.labels}');
        % skipping if channel labels do not match template
        if all(tempChanLabels~=chanNames)
            warning('Check Subject Channel Labels %s %s; skipping\n',lunaid,eeg_date)
            sponpower(subIdx,1) = NaN;
            evopower(subIdx,1) = NaN;
            totpower(subIdx,1) = NaN;     
            continue
        end
        % skip if subject has 0 or 1 epochs (fails if subject has only 1 epoch)
        numEpochs = length(EEG.epoch);
        if numEpochs < 1
            warning('Subject %s %s has %d epochs; skipping\n',lunaid,eeg_date,numEpochs)
            sponpower(subIdx,1) = NaN;
            evopower(subIdx,1) = NaN;
            totpower(subIdx,1) = NaN;            
            continue
        end
        
        % remove baseline from EEG
        EEG = pop_rmbase(EEG,[basestart baseend]);
        
        % define when event starts and ends to reduce size of 4d array for saving
        eventstart_idx = find(EEG.times>=eventstart,1,'first');
        eventend_idx = find(EEG.times>=eventend,1,'first');

        % Compute TF reconstruction of total power
        [sub_totpower, times, freqs] = F_tfpower(EEG,'channels',frow);
        % Just prep period of power
        sub_totpower = sub_totpower(:,eventstart_idx:eventend_idx,:,:);
        % Compute TF reconstruction of evoked power
        [sub_evopower,times,freqs] = F_tfpower(EEG,'channels',frow,'evoked',1);
        sub_evopower = sub_evopower(:,eventstart_idx:eventend_idx,:,:);
        % Compute spontaneous power
        sub_sponpower = sub_totpower - sub_evopower;
        
        % average across all dimensions
        sponpower(subIdx,1) = mean(sub_sponpower,'all');
        totpower(subIdx,1) = mean(sub_totpower,'all');
        evopower(subIdx,1) = mean(sub_evopower,'all');        
    end
    powertable_temp = table(totpower, sponpower, evopower,...
        'VariableNames',...
        {sprintf('%s_broadband_totpower',currtask),...
        sprintf('%s_broadband_sponpower',currtask),...
        sprintf('%s_broadband_evopower',currtask)});
    % join with merged7t_AScore_table
    temptable = [temptable,powertable_temp];
end

%% save table
currdate = datetime();
currdate.Format = 'yyyyMMdd';
broadband_sponevo = temptable;
save(fullfile(resultsdir,['broadband_sponevo_' char(currdate) '.mat']),'broadband_sponevo')
writetable(broadband_sponevo,fullfile(resultsdir,['broadband_sponevo_' char(currdate) '.csv']))