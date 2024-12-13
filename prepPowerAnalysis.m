%% AS Preparatory Period Analysis
%% ERSP and ITC for every channel
% adding necessary paths
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/AS_EEG//PrepPeriodPowerAnalysis/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat') % to get ages

% baseline: [min max] using [-700 -400] from Kai

% Main directory
maindir = hera('Abby/preprocessed_data');
task = 'anti';
taskdirectory = [maindir '/' task]; 
% In path
epochedFolder = [taskdirectory,'/AfterWhole/epochclean_homogenize/'];
epochedName = '*_epochs_kept.set';
% Out path
outFolder = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/TimeFreqPlots/'; % need folder for each participant in this folder

eeglab;

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
    
    % skipping subject if subjectStruct already exisits
    if exist(fullfile(subFolderPath,structName),'file')
        fprintf('skipping; already created %s\n',structName)
        continue
    end

    % loading subject's EEG and defining number of channels
    EEG = pop_loadset(inputfile);
    numChannels = EEG.nbchan;

    % skipping if subject has 0 epochs
    numOfEpochs = length(EEG.epoch);
    if numOfEpochs < 1
        fprintf('skipping; subject has %d epochs\n',numOfEpochs)
        continue
    end

    % getting age from ErrorLatency table
    idxErrorLatTab = find(ErrorLatencyTable.('Luna ID') == str2double(subID) & ErrorLatencyTable.('Scan Date')== str2double(scanDate));
    age = ErrorLatencyTable.Age(idxErrorLatTab);

    % Creating new folder in /TimeFreqPlots/ for each subject to save structure to
    mkdir(outFolder,subFolderName) 
    
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
         subStruct.ERSP.(chanName{j}) = [];
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
    end 

    for chanNum = 1:numChannels
        chanFieldName = sprintf('Ch%d',chanNum);
        
        fprintf('\n--------Processing %s %s Channel %d--------\n',subID,scanDate,chanNum)
        % ersp and itc with no significance level
        % baseline [-700 -400] from Kai
        [ersp, itc, powbase, times,freqs]=newtimef(EEG.data(chanNum,:,:),...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000,EEG.srate,[3 0.5],'freqs',[3 50],...
            'baseline',[-700 -400],'plotersp','off','plotitc','off','plotphasesign','off');
        % saving newtimef() outputs to structures within subStruct for each channel
        subStruct.ERSP.(chanFieldName) = ersp;
        subStruct.ITC.(chanFieldName) = itc;
        subStruct.PowBase.(chanFieldName) = powbase;
        subStruct.Times.(chanFieldName) = times;
        subStruct.Freqs.(chanFieldName) = freqs;

        % ersp and itc plots with significant regions highlighted
        [erspSig, itcSig, powbaseSig,timesSig,freqsSig,erspboot,itcboot]=newtimef(EEG.data(chanNum,:,:),...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000,EEG.srate,[3 0.5],'freqs',[3 50],...
            'baseline',[-700 -400],'alpha',0.05,'plotersp','off','plotitc','off','plotphasesign','off');
        % saving newtimef() outputs to structures within subStruct for each channel
        subStruct.ERSPsig.(chanFieldName) = erspSig;
        subStruct.ITCsig.(chanFieldName) = itcSig;
        subStruct.PowBaseSig.(chanFieldName) = powbaseSig;
        subStruct.TimesSig.(chanFieldName) = timesSig;
        subStruct.FreqsSig.(chanFieldName) = freqsSig;
        subStruct.ERSPbootsig.(chanFieldName) = erspboot;
        subStruct.ITCbootsig.(chanFieldName) = itcboot;
    end
    % save subStruct to subject folder
    save(fullfile(subFolderPath,structName),'subStruct')
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
