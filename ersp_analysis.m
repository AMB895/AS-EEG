% 2025.2.18 amb
%% AS Prep Period Analysis
% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat') % to get ages
% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

%% Setting up directories and paths
% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 

% Uncomment the type of trial you wish to run
trialtype = 1; % correct trials
% trialtype = 0; % incorrect trials
% trialtype = 2; % error corrected trials

% Epoched data path
epochedpath = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];

if trialtype ==1
    epochedname = '*_epochs_kept_cor.set';
    % Save ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials';
    allerspdataname = 'all_ERSP_DATA_cor.mat';
    suberspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/SubjectERSP/';
elseif trialtype ==0
    epochedname = '*_epochs_kept_incor.set';
    % Save ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/IncorrectTrials';
    allerspdataname = 'all_ERSP_DATA_incor.mat';
    suberspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/IncorrectTrials/SubjectERSP/';
elseif trialtype ==2
    epochedname = '*_epochs_kept_errcor.set';
    % Save ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials';
    allerspdataname = 'all_ERSP_DATA_errcor.mat';
    suberspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/SubjectERSP/';
end

% loading template of channel labels
load('/Volumes/Hera/Abby/AS_EEG/templatechannellabels.mat')

%% Computing ERSP for each channel per subject

% Getting all epoched data files
EEGfilenames = dir([epochedpath,epochedname]);
totalEEGs = size(EEGfilenames,1);

% Defining allerspdata and IDmatrix
allerspdata = [];
IDmatrix = [];
% loading allerspdata matrix if it exists
% if it don't exist setting them to empty matricies
% if exist(fullfile(allerspdatapath,allerspdataname),'file')
%     fprintf('loading all ersp data\n')
%     load(fullfile(allerspdatapath,allerspdataname))
% else
%     allerspdata = [];
% end
% 
% % loading ID matrix if it exists
% % if it don't exist setting them to empty matricies
% if exist(fullfile(allerspdatapath,'ID_info.mat'),'file')
%     fprintf('loading ID info\n')
%     load(fullfile(allerspdatapath,'ID_info.mat'))
% else
%     IDmatrix=[];
% end

% loading missing data matrix if exits
% if they don't exist setting them to empty matricies
if exist(fullfile(allerspdatapath,'missingdata.mat'),'file')
    fprintf('loading missing data\n')
    load(fullfile(allerspdatapath,'missingdata.mat'))
else
    notViable=[];
    missingASscore=[];
    noAge=[];
    nochanlabels=[];
    zeroEpochs=[];
end

% setting up progress bar
f = waitbar(0,'Starting','Name','ERSP Progress');

for currentEEG = 1:size(EEGfilenames)
    
    % defining input EEG file from epochclean_homogenize
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedpath,filename];
    [~, currentName, ~ ] = fileparts(inputfile);
    splitCurrentName = split(currentName,'_');
    subID = cell2mat(splitCurrentName(1));
    scanDate = cell2mat(splitCurrentName(2));
    lunaid = [subID,'_',scanDate];
    suberspdataname = sprintf('%sersp.mat',lunaid);
    
    % adding subject to allerspdata matrix if already computed ERSP
    if exist(fullfile(suberspdatapath,suberspdataname),'file')
        fprintf('skipping; computed ersp for %s\nLoading and adding suberspdata to allerspdata\n',lunaid)
        load(fullfile(suberspdatapath,suberspdataname),'suberspdata','times','freqs')
        allerspdata(end+1,:,:,:) = suberspdata;
        fprintf('Adding subject to IDmatrix\n')
        EEG = pop_loadset(inputfile);
        age = EEG.age;
        IDmatrix(end+1,:) = [str2double(subID),str2double(scanDate),age];
        continue
    else

        % defining subject's error latency table
        subIdx = find(ErrorLatencyTable.LunaID==str2double(subID) & ErrorLatencyTable.ScanDate==str2double(scanDate));
        subTable = ErrorLatencyTable(subIdx,:);
        
        % check that subject is viable: has 50% of total trials that are on task
        if subTable.Viable==0
            fprintf('\nSubject %s not viable\n',lunaid)
            notViable(end+1,:) = [str2double(subID),str2double(scanDate)];
            continue
        end
        
        % check that subject has anti score
        if isempty(subTable)
            fprintf('\nSubject %s %s missing from Error Latency Table: no anti score\n',subID,scanDate)
            missingASscore(end+1,:) = [str2double(subID),str2double(scanDate)];
            continue
        end
    
        % loading subject's EEG
        EEG = pop_loadset(inputfile);
        chanNames = string({EEG.chanlocs.labels}');
        numchans = length(chanNames);
    
        % updating progress bar
        waitbar(currentEEG/totalEEGs,f,sprintf('%s %s',subID,scanDate))
    
        % check that subject has an age
        if isfield(EEG,'age')
            age = EEG.age;
        else
            fprintf('\n%s does not have age\n',lunaid)
            noAge(end+1,:) = [str2double(subID),str2double(scanDate)];
            continue
        end
    
        % check that channel names match template channel names
        if all(tempChanLabels~=chanNames)
            fprintf('\nCheck Channel Labels %s\n',lunaid)
            nochanlabels(end+1,:) = [str2double(subID),str2double(scanDate)];
            continue;
        end
    
        % skipping subject if there are zero epochs
        numEpochs = length(EEG.epoch);
        if numEpochs < 1
            fprintf('skipping; %s has %d epochs\n',lunaid,numEpochs)
            zeroEpochs(end+1,:) = [str2double(subID),str2double(scanDate)];
            continue;
        else % calculating ERSP for subjects that have an age, > 0 epochs, are viable, and have an anti score 
    
            % matrix to identify subjects based on index in allerspdata
            IDmatrix(end+1,:) = [str2double(subID),str2double(scanDate),age];
    
            % looping through each channel per subject
            for currentChan = 1:numchans
                % updating progress
                fprintf('\n-------Processing Channel %s Subject %s-------\n',chanNames{currentChan},lunaid)
            
                % Calculating ERSP
                % Inputs to newtimef():
                % Required:
                %   data: 3D array (1,frames,trials)
                %   frames: frames/trial; ignored for 3D data
                %   tlimits: time limits of data epochs, NOT subwindow to extract from epcochs (in ms)
                %   Fs: data sampling rate
                %   varwin: [# of wavelet cycles @ highest freq, window size increase]
                % Optional:
                %   'freqs': frequency limits
                %   'baseline': spectral baseline range [-700 -400] from Kai
                %   'plotersp' & 'plotitc' & 'plotphasesign': plotting options
                % Filtered original data with bandpass filter 0.5-70 Hz -> should not look at frequencies above 70 Hz
                [suberspdata(currentChan,:,:),~,~,times,freqs]=newtimef(EEG.data(currentChan,:,:),...
                     EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate, [3 0.5], 'freqs', [3 70],...
                    'baseline',[-700 -400],'plotersp','off','plotitc','off','plotphasesign','off');   
            end

             % subsetting ersp data for just prep period
            idx580ms = find(times==580);
            suberspdata = suberspdata(:,:,1:idx580ms);
            times = times(1:idx580ms);

            % saving subject's erspdata
            save(fullfile(suberspdatapath,suberspdataname),'suberspdata','times','freqs')

            % add subject's erspdata to all ersp data 4D matrix
            % [subjects,channels,freqs,times]
            allerspdata(end+1,:,:,:) = suberspdata;
            
            clear suberspdata
        end
    end
end
 % Saving ERSP outputs and missing data
save(fullfile(allerspdatapath,allerspdataname),"allerspdata","times","freqs")
save(fullfile(allerspdatapath,'ID_info.mat'),"IDmatrix")
save(fullfile(allerspdatapath,'missingdata.mat'),"zeroEpochs","nochanlabels","noAge","missingASscore","notViable")