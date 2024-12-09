%% AS Preparatory Period Analysis
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
% pop_newtimef() outputs:
% ersp, itc, powbase, times, freqs, erspboot, itcboot, tfdata
%newtimef() inputs:
% itctype: either linear or phase coherence (could be interesting for CFC)
% varwin: [# of cycles, how you want numbers of cycles to increase linearly]
%         default: [3, 0.5] 
% padratio: increases frequency resolution?
%         default: 2
% freqs: [min max] frequency limits
%         default: [minfreq 50] minfreq determined by number of data
%         points, cycles and sampling frequency
% baseline: [min max] using [-700 -400] from Kai
% close all
% clear
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/AS_EEG//PrepPeriodPowerAnalysis/')

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
    filename = [EEGfilenames(currentEEG).name];
    inputfile = [epochedFolder,filename];
    
    [d, currentName, ext ] = fileparts(inputfile);
    splitCurrentName = split(currentName,'_');
    subID = cell2mat(splitCurrentName(1));
    scanDate = cell2mat(splitCurrentName(2));
    subFolderName = [subID,'_',scanDate];
    subFolderPath = [outFolder,subFolderName,'/'];

    mkdir(outFolder,subFolderName) 
    % Need for loop for each channel nested in this for loop, saving the figure for each channel in the participants folder
    EEG = pop_loadset(inputfile);
    numChannels = EEG.nbchan;
    % making files to save newtimef() outputs for all channels
    subStruct = struct('SubID',subID,'ScanDate',scanDate,'ERSP',[], ...
        'ITC',[],'PowBase',[],'Times',[], ...
        'Freqs',[],'ERSPsig',[],'ITCsig',[],...
        'PowBasesig',[],'TimesSig',[],'FreqsSig',[],...
        'ERSPbootsig',[],'ITCbootsig',[]);
    for chanNum = 1:numChannels
        figName = [subID,'_',scanDate,'ChNum',string(chanNum)];
        fig(1)=figure;
        [ersp, itc, powbase, times,freqs]=newtimef(EEG.data(chanNum,:,:),...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000,EEG.srate,[3 0.5],'freqs',[3 50],...
            'baseline',[-700 -400],'plotersp','on','plotitc','on','plotphasesign','off');
        subStruct.ERSP = ersp;
        fig(2)=figure;
        [erspSig, itcSig, powbaseSig,timesSig,freqsSig,erspboot,itcboot]=newtimef(EEG.data(chanNum,:,:),...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000,EEG.srate,[3 0.5],'freqs',[3 50],...
            'baseline',[-700 -400],'alpha',0.05,'plotersp','on','plotitc','on','plotphasesign','off');
        savefig(fig,fullfile(subFolderPath,'TimeFreqPlot.fig'))
    end
end
% get average plots for each channel across adults and children

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
