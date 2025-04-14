%% Calculate ERSP across all subjects/channels and trial types
clear; clc; close all;
% adding necessary paths for all trials
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/AS_EEG/erspFunctions')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

% Setting up directories and paths
% Main directory
maindir = hera('Abby/preprocessed_data');

% Getting age and sex from merge7t
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
merged7t=readtable('merged_7t.csv');

% only care about lunaid, eeg scan date, age at eeg scan and sex from merged7t
merged7t_eeg = [merged7t(:,'lunaid'), merged7t(:,'eeg_date'), merged7t(:,'eeg_age'), merged7t(:,'sex')];

% finding subjects who did not do not have an eeg scan and removing them from merge_7t_eeg
noeeg_idx = find(isnan(merged7t_eeg.eeg_date));
merged7t_eeg(noeeg_idx,:) = [];

merged7t_eeg = renamevars(merged7t_eeg,["lunaid","eeg_date","eeg_age","sex"],["LunaID","ScanDate","Age","Sex"]);

% loading template of channel labels
load('/Volumes/Hera/Abby/AS_EEG/templatechannellabels.mat')

%% Correct AS Trials
% adding path to correct AS trials
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

task = 'anti';
taskdirectory = [maindir '/' task]; 

% Epoched data path and name
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept_cor.set';

% Save ERSP path
corerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials';
corerspdataname = 'corERSPdata.mat';
corsuberspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/SubjectERSP/';
coridmatsavename = 'corIDinfo.csv';

% Getting all epoched data files
EEGfilenames = dir([epochedpath,epochedname]);
totalEEGs = size(EEGfilenames,1);

if exist(fullfile(corerspdatapath,corerspdataname),'file')
    fprintf('Already computed ERSP from all subjects; load from %s/%s\n',corerspdatapath,corerspdataname)
else
    % setting up progress bar
    f = waitbar(0,'Starting','Name','ERSP Progress');
    % Defining allerspdata and IDmatrix and missing data array
    corerspdata = [];
    corIDmatrix = [];
    cormissingdata = [];
    for currentEEG = 1:size(EEGfilenames)
        % defining input into calc_ersp()
        filename = [EEGfilenames(currentEEG).name];
        splitfilename = strsplit(filename,'_');
        subID = cell2mat(splitfilename(1));
        scanDate = cell2mat(splitfilename(2));
        
        % updating waitbar
        waitbar(currentEEG/totalEEGs,f,sprintf('%s %s',subID,scanDate))
        
        % defining subject's error latency table
        subIdx = find(merged7t_eeg.LunaID==str2double(subID) & merged7t_eeg.ScanDate==str2double(scanDate));
        subTable = merged7t_eeg(subIdx,:);
        
        % running calc_ersp for subject
        [corsuberspdata,times,freqs,failcode,subID,scanDate,age,suberspdataname] = calc_ersp(epochedpath,filename,corsuberspdatapath,subTable,tempChanLabels,task);
        if isnan(failcode)
            save(fullfile(corsuberspdatapath,suberspdataname),'corsuberspdata','times','freqs')
            corerspdata(end+1,:,:,:) = corsuberspdata;
            corIDmatrix(end+1,:) = [str2double(subID),str2double(scanDate),age];
        else
            cormissingdata(end+1,:) = [str2double(subID),str2double(scanDate),failcode];
        end
        clear suberspdata
    end
    % Close waitbar
    close(f)
    
    % Save outputs
    save(fullfile(corerspdatapath,corerspdataname),"corerspdata","times","freqs")
    save(fullfile(corerspdatapath,'corIDmatrix.mat'),"corIDmatrix")
    writematrix(corIDmatrix,fullfile(corerspdatapath,coridmatsavename))
    save(fullfile(corerspdatapath,'cormissingdata.mat'),'cormissingdata')
end

%% Error AS Trials
% adding path to correct AS trials
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

task = 'anti';
taskdirectory = [maindir '/' task]; 

% Epoched data path and name
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept_errcor.set';

% Save ERSP path
errcorerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials';
errcorerspdataname = 'errcorERSPdata.mat';
errcorsuberspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/SubjectERSP/';
errcoridmatsavename = 'errcorIDinfo.csv';

% Getting all epoched data files
EEGfilenames = dir([epochedpath,epochedname]);
totalEEGs = size(EEGfilenames,1);

if exist(fullfile(errcorerspdatapath,errcorerspdataname),'file')
    fprintf('Already computed ERSP from all subjects; load from %s/%s\n',errcorerspdatapath,errcorerspdataname)
else
    % setting up progress bar
    f = waitbar(0,'Starting','Name','ERSP Progress');
    
    % Defining allerspdata and IDmatrix and missing data array
    errcorerspdata = [];
    errcorIDmatrix = [];
    errcormissingdata = [];
    for currentEEG = 1:size(EEGfilenames)
        % defining input into calc_ersp()
        filename = [EEGfilenames(currentEEG).name];
        splitfilename = strsplit(filename,'_');
        subID = cell2mat(splitfilename(1));
        scanDate = cell2mat(splitfilename(2));
        
        % updating waitbar
        waitbar(currentEEG/totalEEGs,f,sprintf('%s %s',subID,scanDate))
        
        % defining subject's error latency table
        subIdx = find(merged7t_eeg.LunaID==str2double(subID) & merged7t_eeg.ScanDate==str2double(scanDate));
        subTable = merged7t_eeg(subIdx,:);
        
        % running calc_ersp for subject
        [errcorsuberspdata,times,freqs,failcode,subID,scanDate,age,suberspdataname] = calc_ersp(epochedpath,filename,errcorsuberspdatapath,subTable,tempChanLabels,task);
        if isnan(failcode)
            save(fullfile(errcorsuberspdatapath,suberspdataname),'errcorsuberspdata','times','freqs')
            errcorerspdata(end+1,:,:,:) = errcorsuberspdata;
            errcorIDmatrix(end+1,:) = [str2double(subID),str2double(scanDate),age];
        else
            errcormissingdata(end+1,:) = [str2double(subID),str2double(scanDate),failcode];
        end
        clear suberspdata
    end
    % Close waitbar
    close(f)
    
    % Save outputs
    save(fullfile(errcorerspdatapath,errcorerspdataname),"errcorerspdata","times","freqs")
    save(fullfile(errcorerspdatapath,'errcorIDmatrix.mat'),"errcorIDmatrix")
    writematrix(errcorIDmatrix,fullfile(errcorerspdatapath,errcoridmatsavename))
    save(fullfile(errcorerspdatapath,'errcormissingdata.mat'),'errcormissingdata')
end

%% VGS Trials
% adding path to correct AS trials
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

task = 'vgs';
taskdirectory = [maindir '/' task]; 

% Epoched data path and name
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedname = '*_epochs_kept.set';

% Save ERSP path
vgserspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs';
vgserspdataname = 'vgsERSPdata.mat';
vgssuberspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs/SubjectERSP/';
vgsidmatsavename = 'vgsIDinfo.csv';

% Getting all epoched data files
EEGfilenames = dir([epochedpath,epochedname]);
totalEEGs = size(EEGfilenames,1);

if exist(fullfile(vgserspdatapath,vgserspdataname),'file')
    fprintf('Already computed ERSP from all subjects; load from %s/%s\n',vgserspdatapath,vgserspdataname)
else
    % setting up progress bar
    f = waitbar(0,'Starting','Name','ERSP Progress');
    % Defining allerspdata and IDmatrix and missing data array
    vgserspdata = [];
    vgsIDmatrix = [];
    vgsmissingdata = [];
    for currentEEG = 1:size(EEGfilenames)
        % defining input into calc_ersp()
        filename = [EEGfilenames(currentEEG).name];
        splitfilename = strsplit(filename,'_');
        subID = cell2mat(splitfilename(1));
        scanDate = cell2mat(splitfilename(2));
        
        % updating waitbar
        waitbar(currentEEG/totalEEGs,f,sprintf('%s %s',subID,scanDate))
        
        % defining subject's error latency table
        subIdx = find(merged7t_eeg.LunaID==str2double(subID) & merged7t_eeg.ScanDate==str2double(scanDate));
        subTable = merged7t_eeg(subIdx,:);
        
        % running calc_ersp for subject
        [vgssuberspdata,times,freqs,failcode,subID,scanDate,age,suberspdataname] = calc_ersp(epochedpath,filename,vgssuberspdatapath,subTable,tempChanLabels,task);
        if isnan(failcode)
            save(fullfile(vgssuberspdatapath,suberspdataname),'vgssuberspdata','times','freqs')
            vgserspdata(end+1,:,:,:) = vgssuberspdata;
            vgsIDmatrix(end+1,:) = [str2double(subID),str2double(scanDate),age];
        else
            vgsmissingdata(end+1,:) = [str2double(subID),str2double(scanDate),failcode];
        end
        clear suberspdata
    end
    % Close waitbar
    close(f)
    
    % Save outputs
    save(fullfile(vgserspdatapath,vgserspdataname),"vgserspdata","times","freqs")
    save(fullfile(vgserspdatapath,'vgsIDmatrix.mat'),"vgsIDmatrix")
    writematrix(vgsIDmatrix,fullfile(vgserspdatapath,vgsidmatsavename))
    save(fullfile(vgserspdatapath,'vgsmissingdata.mat'),'vgsmissingdata')
end
