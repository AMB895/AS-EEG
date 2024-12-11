%% AS Preparatory Period Analysis
%% ERSP and ITC for every channel
% adding necessary paths
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/AS_EEG//PrepPeriodPowerAnalysis/')
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat') % to get ages

% baseline: [min max] using [-700 -400] from Kai

% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
epochedFolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedName = '/*_epochs_kept.set';
% Out path
outFolder = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/TimeFreqPlots/'; % need folder for each participant in this folder

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
    name4Struct = sprintf('subject%s%s',subID,scanDate);

    % defining files and folders to be created
    subFolderName = [subID,'_',scanDate];
    subFolderPath = [outFolder,subFolderName,'/'];
    structName = [subID,'_',scanDate,'.mat'];

    % getting age from ErrorLatency table
    idxErrorLatTab = find(ErrorLatencyTable.('Luna ID') == subID & ErrorLatencyTable.('Scan Date')==scanDate);
    age = ErrorLatencyTable(idxErrorLatTab,'Age');

    % Creating new folder in /TimeFreqPlots/ for each subject to save channel plots to
    mkdir(outFolder,subFolderName) 

    % loading subject's EEG and defining number of channels
    EEG = pop_loadset(inputfile);
    numChannels = EEG.nbchan;
    
    % making structure to save newtimef() outputs for each channel
    subStruct = struct('SubID',subID,'ScanDate',scanDate,'Age',age,'ERSP',[], ...
        'ITC',[],'PowBase',[],'Times',[], ...
        'Freqs',[],'ERSPsig',[],'ITCsig',[],...
        'PowBaseSig',[],'TimesSig',[],'FreqsSig',[],...
        'ERSPbootsig',[],'ITCbootsig',[]);
    
    % naming channels for structures within subStruct
    for i=1:numChannels
        chanName{i}= sprintf('Ch%d',i);
    end

    % preallocating for number of channels to each field in subStruct
    for j = 1:length(chanName)
         subStruct.ERSP.(chanName{j}) = []; % Dynamically create a field
         subStruct.ITC.(chanName{j}) = [];
         subStruct.PowBase.(chanName{j}) = [];
         subStruct.Times.(chanName{j}) = [];
         subStruct.Freqs.(chanName{j}) = [];
         subStruct.ERSPsig.(chanName{j}) = [];
         subStruct.ITCsig.(chanName{j}) = [];
         subStruct.PowBaseSig.(chanName{j}) = [];
         subStruct.TimesSig.(chanName{j}) = [];
         subStruct.FreqsSig.(chanName{j}) = [];
         subStruct.ERSPbootsig.(chanName{j}) = [];
         subStruct.ITCbootsig.(chanName{j}) = [];
    % preallocating strucutres to save all ersp, itc, times and freqs
         if subStruct.Age >= 18
            allerspAdult.(name4Struct).(chanName{j}) = [];
            allitcAdult.(name4Struct).(chanName{j}) = [];
            alltimesAdult.(name4Struct).(chanName{j}) = [];
            allfreqsAdult.(name4Struct).(chanName{j}) = [];
            allerspsigAdult.(name4Struct).(chanName{j}) = [];
            allitcsigAdult.(name4Struct).(chanName{j}) = [];
            alltimessigAdult.(name4Struct).(chanName{j}) = [];
            allfreqssigAdult.(name4Struct).(chanName{j}) = [];
            allerspbootAdult.(name4Struct).(chanName{j}) = [];
            allitcbootAdult.(name4Struct).(chanName{j}) = [];
         elseif subStruct.Age < 18
            allerspChild.(name4Struct).(chanName{j}) = [];
            allitcChild.(name4Struct).(chanName{j}) = [];
            alltimesChild.(name4Struct).(chanName{j}) = [];
            allfreqsChild.(name4Struct).(chanName{j}) = [];
            allerspsigChild.(name4Struct).(chanName{j}) = [];
            allitcsigChild.(name4Struct).(chanName{j}) = [];
            allerspbootChild.(name4Struct).(chanName{j}) = [];
            allitcbootChild.(name4Struct).(chanName{j}) = [];
         end
    end 

    % Checking if subject structure file is already created
    if exist(fullfile(subFolderPath,structName),'file')
        fprintf('skipping; already created %s\n',structName)
        continue
    end

    for chanNum = 1:numChannels
        % naming conventions for saving figures and outputs from newtimef()
        figName = [subID,'_',scanDate,'_Ch',string(chanNum)];
        chanFieldName = sprintf('Ch%d',chanNum);

        % Checking to see if figure has already been created for channel
        figFile = [fullfile(subFolderPath,figName),'.fig'];
        if exist(figFile,'file' ) % if this doesn't work maybe do not use 'file'
            fprintf('skipping; already created %s\n',figName)
            continue
        end

        fig(1)=figure;  % fig contains both figures created for each channel
        % ersp and itc plots with no significance level
        % baseline [-700 -400] from Kai
       
        % outputs from newtimef() MUST be the same length for each subject to get averages
        [ersp, itc, powbase, times,freqs]=newtimef(EEG.data(chanNum,:,:),...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000,EEG.srate,[3 0.5],'freqs',[3 50],...
            'baseline',[-700 -400],'plotersp','on','plotitc','on','plotphasesign','off');
        % saving newtimef() outputs to structures within subStruct for each channel
        subStruct.ERSP.(chanFieldName) = ersp;
        subStruct.ITC.(chanFieldName) = itc;
        subStruct.PowBase.(chanFieldName) = powbase;
        subStruct.Times.(chanFieldName) = times;
        subStruct.Freqs.(chanFieldName) = freqs;

        fig(2)=figure;
        % ersp and itc plots with significant regions highlighted
        % outputs from newtimef() MUST be the same length for each subject to get averages
        [erspSig, itcSig, powbaseSig,timesSig,freqsSig,erspboot,itcboot]=newtimef(EEG.data(chanNum,:,:),...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000,EEG.srate,[3 0.5],'freqs',[3 50],...
            'baseline',[-700 -400],'alpha',0.05,'plotersp','on','plotitc','on','plotphasesign','off');
        % saving newtimef() outputs to structures within subStruct for each channel
        subStruct.ERSPsig.(chanFieldName) = erspSig;
        subStruct.ITCsig.(chanFieldName) = itcSig;
        subStruct.PowBaseSig.(chanFieldName) = powbaseSig;
        subStruct.TimesSig.(chanFieldName) = timesSig;
        subStruct.FreqsSig.(chanFieldName) = freqsSig;
        subStruct.ERSPbootsig.(chanFieldName) = erspboot;
        subStruct.ITCbootsig.(chanFieldName) = itcboot;
        % saving outputs into specific age structure
        if subStruct.Age >= 18
            allerspAdult.(name4Struct).(chanFieldName) = ersp;
            allitcAdult.(name4Struct).(chanFieldName) = itc;
            alltimesAdult.(name4Struct).(chanFieldName) = times;
            allfreqsAdult.(name4Struct).(chanFieldName) = freqs;
            allerspsigAdult.(name4Struct).(chanFieldName) = erspSig;
            allitcsigAdult.(name4Struct).(chanFieldName) = itcSig;
            alltimessigAdult.(name4Struct).(chanFieldName) = timesSig;
            allfreqssigAdult.(name4Struct).(chanFieldName) = freqsSig;
            allerspbootAdult.(name4Struct).(chanFieldName) = erspboot;
            allitcbootAdult.(name4Struct).(chanFieldName) = itcboot;
        elseif subStruct.Age < 18
            allerspChild.(name4Struct).(chanFieldName) = ersp;
            allitcChild.(name4Struct).(chanFieldName) = itc;
            alltimesChild.(name4Struct).(chanFieldName) = times;
            allfreqsChild.(name4Struct).(chanFieldName) = freqs;
            allerspsigChild.(name4Struct).(chanFieldName) = erspSig;
            allitcsigChild.(name4Struct).(chanFieldName) = itcSig;
            alltimessigChild.(name4Struct).(chanFieldName) = timesSig;
            allfreqssigChild.(name4Struct).(chanFieldName) = freqsSig;
            allerspbootChild.(name4Struct).(chanFieldName) = erspboot;
            allitcbootChild.(name4Struct).(chanFieldName) = itcboot;
        end
        % saving both figures in fig in the subject folder
        savefig(fig,fullfile(subFolderPath,figName)) 
        close all
        clear fig
    end
    % save subStruct to subject folder
    save(fullfile(subFolderPath,structName),'subStruct')
end

% Path for all data structures
adultDataPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Adults/';
childDataPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Children/';
allDataAdult = 'TimeFreqDataAdult.mat';
allDataChildren = 'TimeFreqDataChildren.mat';
% save all data structures to one .mat file
if subStruct.Age >= 18    
    save(fullfile(adultDataPath,allDataAdult),'allerspAdult','allitcAdult','alltimesAdult',...
        'allfreqsAdult','allerspsigAdult','allitcsigAdult','alltimessigAdult','allfreqssigAdult',...
        'allerspbootAdult','allitcbootAdult')
elseif subStruct.Age < 18
    save(fullfile(childDataPath,allDataChildren),'allerspChild','allitcChild','alltimesChild',...
        'allfreqsChild','allerspsigChild','allitcsigChild','alltimessigChild','allfreqssigChild',...
        'allerspbootChild','allitcbootChild')
end
%% Averaging ERSP and ITC across children and adults for all channels
% load allersp, allitc, etc data structures for children and adults
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Adults/TimeFreqDataAdult.mat'); 
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Children/TimeFreqDataChildren.mat');
numAdults = length(fieldnames(allerspAdult)); % number of adults
adultFieldNames = fieldnames(allerspAdult); % names of adults
numChild = length(fieldnames(allerspChild)); % number of children
childFieldNames = fieldnames(allerspChild); % names of children

% Paths to save channel data to
adultPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Adults/';
childPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Children/';

% getting field names for channels
numOfChan = 64; %is there actually 64 channels for every subject?
for i=1:numOfChan
    chanName{i}= sprintf('Ch%d',i);
end

%% averaging channel data from adults
for h = 1:numOfChan
    channelfieldname = chanName{h};
    % defining file names to name averages for each channel
    erspFileName = sprintf('adultERSP%s.mat',channelfieldname);
    itcFileName = sprintf('adultITC%s.mat',channelfieldname);
    timesFileName = sprintf('adultTimes%s.mat',channelfieldname);
    freqsFileName = sprintf('adultFreqs%s.mat',channelfieldname);
    erspsigFileName = sprintf('adultERSPsig%s.mat',channelfieldname);
    itcsigFileName = sprintf('adultITCsig%s.mat',channelfieldname);
    timessigFileName = sprintf('adultTimessig%s.mat',channelfieldname);
    freqssigFileName = sprintf('adultFreqssig%s.mat',channelfieldname);
    erspbootFileName = sprintf('adultERSPboot%s.mat',channelfieldname);
    itcbootFileName = sprintf('adultITCboot%s.mat',channelfieldname);
    for n = 1:numAdults
        % looping through each subject to get data and storing each channel's data in new variable
        subjectfieldname = adultFieldNames{n};
        erspChanData = allerspAdult.(subjectfieldname).(channelfieldname);
        erspAllSubChanAdult(:,:,n) = erspChanData; % 3rd dimension is subject
        itcChanData = allitcAdult.(subjectfieldname).(channelfieldname);
        itcAllSubChanAdult(:,:,n) = itcChanData;
        timesChanData = alltimesAdult.(subjectfieldname).(channelfieldname);
        timesAllSubChanAdult(:,:,n) = timesChanData;
        freqsChanData = allfreqsAdult.(subjectfieldname).(channelfieldname);
        freqsAllSubChanAdult(:,:,n) = freqsChanData;

        % data from setting significance level
        erspsigChanData = allerspsigAdult.(subjectfieldname).(channelfieldname);
        erspsigAllSubChanAdult(:,:,n) = erspsigChanData;
        itcsigChanData = allitcsigAdult.(subjectfieldname).(channelfieldname);
        itcsigAllSubChanAdult(:,:,n) = itcsigChanData;
        timessigChanData = alltimessigAdult.(subjectfieldname).(channelfieldname);
        timessigAllSubChanAdult(:,:,n) = timessigChanData;
        freqssigChanData = allfreqssigAdult.(subjectfieldname).(channelfieldname);
        freqssigAllSubChanAdult(:,:,n) = freqssigChanData;
        erspbootChanData = allerspbootAdult.(subjectfieldname).(channelfieldname);
        erspbootAllSubChanAdult(:,:,n) = erspbootChanData; % 3rd dimension is subject
        itcbootChanData = allitcsigAdult.(subjectfieldname).(channelfieldname);
        itcbootAllSubChanAdult(:,:,n) = itcbootChanData;
    end
    % before looping to next channel get averages for ersp and itc (maybe
    % times and freqs but they should be the same for each channel and
    % subject)
    avgERSP_adults = mean(erspAllSubChanAdult,3);
    avgITC_adults = mean(itcAllSubChanAdult,3);
    avgTimes_adults = mean(timesAllSubChanAdult,3);
    avgFreqs_adults = mean(freqsAllSubChanAdult,3);
    avgERSPsig_adults = mean(erspsigAllSubChanAdult,3);
    avgITCsig_adults = mean(itcsigAllSubChanAdult,3);
    avgTimessig_adults = mean(timessigAllSubChanAdult,3);
    avgFreqssig_adults = mean(freqssigAllSubChanAdult,3);
    avgERSPboot_adults = mean(erspbootAllSubChanAdult,3);
    avgITCboot_adults = mean(itcbootAllSubChanAdult,3);

    % Plotting and saving plots for average of every channel
    % ERSP and ITC average plot
    Fig = figure;
    subplot(2,1,1)
    tftopo(avgERSP_adults,avgTimes_adults,avgFreqs_adults)
    subplot(2,1,2)
    tftopo(avgITC_adults,avgTimes_adults,avgFreqs_adults)

    figName = sprintf('AverageERSP/ITC_Adults%s',channelfieldname);
    savefig(Fig,fullfile(adultPath,figName))
    close all
    clear Fig

    % not sure exactly how to do significance levels plots- need to mess
    % around with that

    % saving average across subjects for each channel
    save(fullfile(adultPath,erspFileName),'avgERSP_adults')
    save(fullfile(adultPath,itcFileName),'avgITC_adults')
    save(fullfile(adultPath,timesFileName),'avgTimes_adults')
    save(fullfile(adultPath,freqsFileName),'avgFreqs_adults')
    save(fullfile(adultPath,erspsigFileName),'avgERSPsig_adults')
    save(fullfile(adultPath,itcsigFileName),'avgITCsig_adults')
    save(fullfile(adultPath,timessigFileName),'avgTimessig_adults')
    save(fullfile(adultPath,freqssigFileName),'avgFreqssig_adults')
    save(fullfile(adultPath,erspbootFileName),'avgERSPboot_adults')
    save(fullfile(adultPath,itcbootFileName),'avgITCboot_adults')
    % maybe instead of saving every subjects channel data, average all
    % subjects once the loop goes through all the subjects
    % save(fullfile(adultPath,erspFileName),"erspAllSubChanAdult")
    % save(fullfile(adultPath,itcFileName),"itcAllSubChanAdult")
    % save(fullfile(adultPath,timesFileName),"timesAllSubChanAdult")
    % save(fullfile(adultPath,freqsFileName),"freqsAllSubChanAdult")
end

%% averaging channel data for children
for h = 1:numOfChan
    channelfieldname = chanName{h};
    % defining file names to name averages for each channel
    erspFileName = sprintf('childERSP%s.mat',channelfieldname);
    itcFileName = sprintf('childITC%s.mat',channelfieldname);
    timesFileName = sprintf('childTimes%s.mat',channelfieldname);
    freqsFileName = sprintf('childFreqs%s.mat',channelfieldname);
    erspsigFileName = sprintf('childERSPsig%s.mat',channelfieldname);
    itcsigFileName = sprintf('childITCsig%s.mat',channelfieldname);
    timessigFileName = sprintf('childTimessig%s.mat',channelfieldname);
    freqssigFileName = sprintf('childFreqssig%s.mat',channelfieldname);
    erspbootFileName = sprintf('childERSPboot%s.mat',channelfieldname);
    itcbootFileName = sprintf('childITCboot%s.mat',channelfieldname);
    for n = 1:numChild
         % looping through each subject to get data and storing each channel's data in new variable
        subjectfieldname = childFieldNames{n};
        erspChanData = allerspChild.(subjectfieldname).(channelfieldname);
        erspAllSubChanChild(:,:,n) = erspChanData; % 3rd dimension is subject
        itcChanData = allitcChild.(subjectfieldname).(channelfieldname);
        itcAllSubChanChild(:,:,n) = itcChanData;
        timesChanData = alltimesChild.(subjectfieldname).(channelfieldname);
        timesAllSubChanChild(:,:,n) = timesChanData;
        freqsChanData = allfreqsChild.(subjectfieldname).(channelfieldname);
        freqsAllSubChanChild(:,:,n) = freqsChanData;

        % data from setting significance level
        erspsigChanData = allerspsigChild.(subjectfieldname).(channelfieldname);
        erspsigAllSubChanChild(:,:,n) = erspsigChanData;
        itcsigChanData = allitcsigChild.(subjectfieldname).(channelfieldname);
        itcsigAllSubChanChild(:,:,n) = itcsigChanData;
        timessigChanData = alltimessigChild.(subjectfieldname).(channelfieldname);
        timessigAllSubChanChild(:,:,n) = timessigChanData;
        freqssigChanData = allfreqssigChild.(subjectfieldname).(channelfieldname);
        freqssigAllSubChanChild(:,:,n) = freqssigChanData;
        erspbootChanData = allerspbootChild.(subjectfieldname).(channelfieldname);
        erspbootAllSubChanChild(:,:,n) = erspbootChanData; % 3rd dimension is subject
        itcbootChanData = allitcsigChild.(subjectfieldname).(channelfieldname);
        itcbootAllSubChanChild(:,:,n) = itcbootChanData;
    end
    % before looping to next channel get averages for ersp and itc (maybe
    % times and freqs but they should be the same for each channel and
    % subject)
    avgERSP_child = mean(erspAllSubChanChild,3);
    avgITC_child = mean(itcAllSubChanChild,3);
    avgTimes_child = mean(timesAllSubChanChild,3);
    avgFreqs_child = mean(freqsAllSubChanChild,3);
    avgERSPsig_child = mean(erspsigAllSubChanChild,3);
    avgITCsig_child = mean(itcsigAllSubChanChild,3);
    avgTimessig_child = mean(timessigAllSubChanChild,3);
    avgFreqssig_child = mean(freqssigAllSubChanChild,3);
    avgERSPboot_child = mean(erspbootAllSubChanChild,3);
    avgITCboot_child = mean(itcbootAllSubChanChild,3);

    % Plotting and saving plots for average of every channel
    % ERSP and ITC average plot
    Fig = figure;
    subplot(2,1,1)
    tftopo(avgERSP_child,avgTimes_child,avgFreqs_child)
    subplot(2,1,2)
    tftopo(avgITC_child,avgTimes_child,avgFreqs_child)

    figName = sprintf('AverageERSP/ITC_Children%s',channelfieldname);
    savefig(Fig,fullfile(childPath,figName))
    close all
    clear Fig

    % not sure exactly how to do significance levels plots- need to mess
    % around with that

    % saving average across subjects for each channel
    save(fullfile(childPath,erspFileName),'avgERSP_child')
    save(fullfile(childPath,itcFileName),'avgITC_child')
    save(fullfile(childPath,timesFileName),'avgTimes_child')
    save(fullfile(childPath,freqsFileName),'avgFreqs_child')
    save(fullfile(childPath,erspsigFileName),'avgERSPsig_child')
    save(fullfile(childPath,itcsigFileName),'avgITCsig_child')
    save(fullfile(childPath,timessigFileName),'avgTimessig_child')
    save(fullfile(childPath,freqssigFileName),'avgFreqssig_child')
    save(fullfile(childPath,erspbootFileName),'avgERSPboot_child')
    save(fullfile(childPath,itcbootFileName),'avgITCboot_child')
    % maybe instead of saving every subjects channel data, average all
    % subjects once the loop goes through all the subjects
    % save(fullfile(adultPath,erspFileName),"erspAllSubChanAdult")
    % save(fullfile(adultPath,itcFileName),"itcAllSubChanAdult")
    % save(fullfile(adultPath,timesFileName),"timesAllSubChanAdult")
    % save(fullfile(adultPath,freqsFileName),"freqsAllSubChanAdult")
end


%% used to get % change in baseline from signal
% this will be used later maybe
% Baseline from epochs to compare to prep period
% Baseline 700 to 400 ms before cue presentation (Hwang 2016)
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Baseline/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/PrepPeriod/')
% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
epochedFolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedName = '/*_epochs_rj.set';
% Out path
baselineFilePath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Baseline';
prepFilePath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/PrepPeriod';

EEGfilenames = dir([epochedFolder,epochedName]);
for currentEEG=1:size(EEGfilenames)
    % Baseline epoch 700 ms to 400 ms before red fixation cross
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedFolder,filename];
    [d, currentName, ext ] = fileparts(inputfile);
    baselineName = [currentName,'_epoch_baseline'];
    baselineFile = [baselineName,'.set'];
    if exist(baselineFile,'file')
        fprintf('skipping; already created %s\n',baselineName)
    else
        EEG = pop_loadset(inputfile);
        fprintf('Baseline epoching %s\n',baselineName)
        baselineEEG = pop_select(EEG,'time',[-0.7 -0.4]);
        baselineEEG = pop_saveset(baselineEEG,'filename',baselineName,'filepath',baselineFilePath);
    end

    % Smaller prep period epoch -100 ms before red fixation cross to 500 ms after red fixation cross
    prepName = [currentName,'_epoch_prep'];
    prepFile = [prepName,'.set'];
    if exist(prepFile,'file')
        fprintf('skipping; already created %s\n',prepName)
    else
        EEG = pop_loadset(inputfile);
        fprintf('Prep Period epoching %s\n',prepName)
        prepEEG = pop_select(EEG,'time',[-0.1 0.5]);
        prepEEG = pop_saveset(prepEEG,'filename',prepName,'filepath',prepFilePath);
    end
end
