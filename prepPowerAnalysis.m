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
        % may need to specify points and time window bc equal size matricies will be needed to calculate averages across individuals
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
        % ersp and itc plots with significant regions
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
%% Creating average ERSP and ITC plots for Children and Adults for all channels
adultData = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Adults/TimeFreqDataAdult.mat'); % loads allersp, etc data structures
childData = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Children/TimeFreqDataChildren.mat');
numAdults = length(fieldnames(allerspAdult)); % number of adults
adultFieldNames = fieldnames(allerspAdult); % names of adults
numChild = length(fieldnames(allerspChild)); % number of children
childFieldNames = fieldnames(allerspChild); % names of children

% Paths to save channel data to
adultPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Adults/';
childPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Children/';

% getting field names of channels
numOfChan = 64; %is there actually 64 channels for every subject?
for i=1:numOfChan
    chanName{i}= sprintf('Ch%d',i);
end

% getting channel data from adults
for h = 1:numOfChan
    channelfieldname = chanName{h};
    erspFileName = sprintf('adultERSP%s.mat',channelfieldname);
    itcFileName = sprintf('adultITC%s.mat',channelfieldname);
    timesFileName = sprintf('adultTimes%s.mat',channelfieldname);
    freqsFileName = sprintf('adultFreqs%s.mat',channelfieldname);
    for n = 1:numAdults
        subjectfieldname = adultFieldNames{n};
        erspChanData = allerspAdult.(subjectfieldname).(channelfieldname);
        erspAllSubChanAdult(:,:,n) = erspChanData; % 3rd dimension is subject
        itcChanData = allitcAdult.(subjectfieldname).(channelfieldname);
        itcAllSubChanAdult(:,:,n) = itcChanData;
        timesChanData = alltimesAdult.(subjectfieldname).(channelfieldname);
        timesAllSubChanAdult(:,:,n) = timesChanData;
        freqsChanData = allfreqsAdult.(subjectfieldname).(channelfieldname);
        freqsAllSubChanAdult(:,:,n) = freqsChanData;
    end
    save(fullfile(adultPath,erspFileName),"erspAllSubChanAdult")
    save(fullfile(adultPath,itcFileName),"itcAllSubChanAdult")
    save(fullfile(adultPath,timesFileName),"timesAllSubChanAdult")
    save(fullfile(adultPath,freqsFileName),"freqsAllSubChanAdult")
end

% outputs from newtimef() MUST be the same length for each subject to get averages


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
