%% Calculate time/frequency reconstruction and ERSP across all subjects/channels and trial types
%% On PSC
clear; clc; close all;
if exist('/Volumes/Hera','dir')
    error('Run on PSC')
elseif exist('/ocean/projects/','dir')
    fprintf('Running Time-Frequency ERSP Analysis on PSC\n')
    
    % add paths to tool boxes and main data directory
    addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/tools/')
    addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/')
    
    % add path to function that calculates ERSP, spon, and evoked power
    addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/functions/')
    
    % start eeglab
    [ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;
    
    % Names and paths of data to be loaded
    merge7t_eeg = readtable('/ocean/projects/soc230004p/shared/antisaccade_eeg/merge7t_eeg.csv');

    % finding subjects who did not do not have an eeg scan and removing them from merge7t_eeg
    noeeg_idx = find(isnan(merge7t_eeg.eeg_date));
    merge7t_eeg(noeeg_idx,:) = [];

    % loading template of channel labels
    load('/ocean/projects/soc230004p/shared/antisaccade_eeg/templatechannellabels.mat')

    % adding path to correct AS trials and VGS trials for prep period
    preprocesseddir_prep='/ocean/projects/soc230004p/shared/antisaccade_eeg/preprocessed_data/prep/';
    preprocesseddir_stim = '/ocean/projects/soc230004p/shared/antisaccade_eeg/preprocessed_data/stim/';    
    % FIX ME!! add stim onset epochs and fix folders and rsync stimonset eegs
    
    % Add all eeg files to one big structure
    EEGfilenames = {dir([preprocesseddir_prep,'*_epochs_kept_cor.set']),...
        dir([preprocesseddir_prep,'*_epochs_kept.set']),...
        dir([preprocesseddir_prep,'*_epochs_kept_errcor.set']),...
        dir([preprocesseddir_stim,'*_epochs_kept_cor.set']),...
        dir([preprocesseddir_stim,'*_epochs_kept.set'])};
    tasks_str = ["CorAS","VGS","ErrCorAS","CorASstim","VGSstim"];

for currentTask = 1:length(tasks_str)
    % defining current task and total eegs for that task
    taskfilenames = EEGfilenames{currentTask};
    tasktotalEEGs = length(taskfilenames);
    currenttask_str = tasks_str(currentTask);
    
    % all subjects data save path for current task
    tasksavepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/%s/',currenttask_str);
    ersptasksavename = sprintf('%sersp.mat',currenttask_str);
    erstasksavename = sprintf('%sers.mat',currenttask_str);
    evokedtasksavename = sprintf('%sevokedpower.mat',currenttask_str);
    spontaneoustasksavename = sprintf('%sspontaneouspower.mat',currenttask_str);
    idmatrixtasksavename = sprintf('%sidmatrix.mat',currenttask_str);
    missingdatatasksavename = sprintf('%smissingdata.mat',currenttask_str);
    
    % determine if task's full data set was saved
    if exist(fullfile(tasksavepath,ersptasksavename),'file')
        fprintf('computed ERSP and time/freq data for %s\n',currenttask_str)
    else
        allersp = [];
        allevokedpreppower = [];
        allspontaneouspreppower = [];
        allers = [];
        missingdata = [];
        idmatrix = [];
        for currentEEG = 1:tasktotalEEGs
            % getting subject's id and scan date
            filename = [taskfilenames(currentEEG).name];
            filepath = [taskfilenames(currentEEG).folder];
            splitfilename = strsplit(filename,'_');
            subID = cell2mat(splitfilename(1));
            scanDate = cell2mat(splitfilename(2));
            
            % subject data save path
            subtasksavepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/%s/SubjectData/',currenttask_str);
            subtasksavename = sprintf('%s_%s_%sersp.mat',subID,scanDate,currenttask_str);
            
            % determine if subject's ersp and time/freq has already been computed
            if exist(fullfile(subtasksavepath,subtasksavename),'file')
                fprintf('skipping; computed ERSP and time/freq for %s_%s; loading file\n',subID,scanDate)
                load(fullfile(subtasksavepath,subtasksavename))
                allersp(end+1,:,:,:) = sub_ersp;
                allevokedpreppower(end+1,:,:,:) = sub_evokedpreppower;
                allspontaneouspreppower(end+1,:,:,:) = sub_spontaneouspreppower;
                allers(end+1,:,:,:) = sub_ers;
                % defining subject's error latency table to get age
                subIdx = find(merge7t_eeg.lunaid==str2double(subID) & merge7t_eeg.eeg_date==str2double(scanDate));
                age = table2array(merge7t_eeg(subIdx,'eeg_age'));
                % adding subject to idmatrix
                idmatrix(end+1,:) = [str2double(subID) str2double(scanDate) age];
            else
                % defining subject's error latency table
                subIdx = find(merge7t_eeg.lunaid==str2double(subID) & merge7t_eeg.eeg_date==str2double(scanDate));
                subTable = merge7t_eeg(subIdx,:);
                
                % skipping if not in Merge 7t
                if isempty(subTable)
                    fprintf('Subject %s %s is not in merge 7t; skipping\n',subID,scanDate)
                    failcode = 1;
                    missingdata(end+1,:) = [str2double(subID) str2double(scanDate) failcode];
                    continue
                end
                
                % load subject's EEG
                EEG = pop_loadset('filename',filename,'filepath',filepath);
                chanNames = string({EEG.chanlocs.labels}');

                % getting age from EEG structure
                if isfield(EEG,'age')
                    age = EEG.age;
                else % skipping if subject does not have age
                    fprintf('Subject %s %s does not have age; skipping\n',subID,scanDate)
                    failcode = 2;
                    missingdata(end+1,:) = [str2double(subID) str2double(scanDate) failcode];
                    continue;
                end
                
                % skipping if channel labels do not match template
                % check that channel labels are the same
                if all(tempChanLabels~=chanNames)
                    fprintf('Check Subject Channel Labels %s %s; skipping\n',subID,scanDate)
                    failcode = 3;
                    missingdata(end+1,:) = [str2double(subID) str2double(scanDate) failcode];
                    continue;
                end
                
                % skip if subject has zero epochs
                numEpochs = length(EEG.epoch);
                if numEpochs < 1
                    fprintf('Subject %s %s has %d epochs; skipping\n',subID,scanDate,numEpochs)
                    failcode = 0;
                    missingdata(end+1,:) =  [str2double(subID) str2double(scanDate) failcode];
                    continue;
                end
                
                % run calc_ersp_timefreq()
                [sub_ersp,sub_ers,sub_evokedpreppower,sub_spontaneouspreppower,times,freqs,preptimes] = calc_ersp_timefreq(EEG);
                idmatrix(end+1,:) = [str2double(subID) str2double(scanDate) age];
                
                % save subject's data
                save(fullfile(subtasksavepath,subtasksavename),'sub_ersp','sub_ers','sub_evokedpreppower',...
                    'sub_spontaneouspreppower','times','freqs','preptimes')
                
                % adding subject's data to full matrix
                allersp(end+1,:,:,:) = sub_ersp;
                allers(end+1,:,:,:) = sub_ers;
                allevokedpreppower(end+1,:,:,:) = sub_evokedpreppower;
                allspontaneouspreppower(end+1,:,:,:) = sub_spontaneouspreppower;
           end
            clear sub_ersp sub_ers sub_evokedpreppower sub_spontaneouspreppower
        end
        
        % save ersp, ers, spontaneous and evoked in separate files with prep times and freqs
        save(fullfile(tasksavepath,ersptasksavename),'allersp','preptimes','freqs')
        save(fullfile(tasksavepath,erstasksavename),'allers','times','freqs')
        save(fullfile(tasksavepath,evokedtasksavename),'allevokedpreppower','preptimes','freqs')
        save(fullfile(tasksavepath,spontaneoustasksavename),'allspontaneouspreppower','preptimes','freqs')
        save(fullfile(tasksavepath,idmatrixtasksavename),'idmatrix')
        save(fullfile(tasksavepath,missingdatatasksavename),'missingdata')
        
        % clear data
        clear allersp allers allevokedpreppower allspontaneouspreppower preptimes times freqs idmatrix missingdata 
    end     
end
end 