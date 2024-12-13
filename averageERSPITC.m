%% getting average plots for adults and children
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat') % to get ages
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/functions/miscfunc')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath(genpath('/resources/Euge/'))
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/')
% In path
allTimeFreqFolder = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/TimeFreqPlots/';
allSubFolders = dir(allTimeFreqFolder);
allSubFolders(2) = []; % first two entries are blank?
allSubFolders(1) = [];
% Outpaths
adultDataPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Adults/';
childDataPath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/Children/';

% Splitting up children and adults
for chan=1:64 
    chanFieldName= sprintf('Ch%d',chan);
    % names of files being saved
    adultAvgName = sprintf('AdultAverages%s.mat',chanFieldName);
    childAvgName = sprintf('ChildAverages%s.mat',chanFieldName);
    % adultERSPname = sprintf('adultERSP%s',chanFieldName);
    % adultITCname = sprintf('adultITC%s',chanFieldName);
    % adultTimesname = sprintf('adultTimes%s',chanFieldName);
    % adultFreqsname = sprintf('adultFreqs%s',chanFieldName);
    % adultPowBasename = sprintf('adultPowBase%s',chanFieldName);
    % adultERSPsigname = sprintf('adultERSPsig%s',chanFieldName);
    % adultITCsigname = sprintf('adultITCsig%s',chanFieldName);
    % adultPowBasesigname = sprintf('adultPowBasesig%s',chanFieldName);
    % childERSPname = sprintf('childERSP%s',chanFieldName);
    % childITCname = sprintf('childITC%s',chanFieldName);
    % childTimesname = sprintf('childTimes%s',chanFieldName);
    % childFreqsname = sprintf('childFreqs%s',chanFieldName);
    % childPowBasename = sprintf('childPowBase%s',chanFieldName);
    % childERSPsigname = sprintf('childERSPsig%s',chanFieldName);
    % childITCsigname = sprintf('childITCsig%s',chanFieldName);
    % childPowBasesigname = sprintf('childPowBasesig%s',chanFieldName);
    a = 1;
    c=1;
    for currentSub = 1:size(allSubFolders)
        subjectFolderName = sprintf('%s',allSubFolders(currentSub).name);
        structureName = sprintf('%s.mat',allSubFolders(currentSub).name);
        subFieldName = sprintf('sub%s',subjectFolderName);
        % loaded variable is called subStruct
        load(sprintf('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/TimeFreqPlots/%s/%s',subjectFolderName,structureName));
        age = subStruct.Age;
        if age >= 18
            % do these have to be cell arrays?
            adultERSP(:,:,a)=subStruct.ERSP.(chanFieldName); % 3rd dimension is subject
            % newa = cell2mat(adultERSP);
            adultITC(:,:,a)=subStruct.ITC.(chanFieldName);
            adultTimes(:,:,a)=subStruct.Times.(chanFieldName);
            adultFreqs(:,:,a)=subStruct.Freqs.(chanFieldName);
            adultPowBase(:,:,a)=subStruct.PowBase.(chanFieldName);
            adultERSPsig(:,:,a)=subStruct.ERSPsig.(chanFieldName);
            adultITCsig(:,:,a)=subStruct.ITCsig.(chanFieldName);
            adultPowBasesig(:,:,a)=subStruct.PowBaseSig.(chanFieldName);
            a=a+1;
        elseif age < 18
            childERSP(:,:,c)=subStruct.ERSP.(chanFieldName);
            childITC(:,:,c)=subStruct.ITC.(chanFieldName);
            childTimes(:,:,c)=subStruct.Times.(chanFieldName);
            childFreqs(:,:,c)=subStruct.Freqs.(chanFieldName);
            childPowBase(:,:,c)=subStruct.PowBase.(chanFieldName);
            childERSPsig(:,:,c)=subStruct.ERSPsig.(chanFieldName);
            childITCsig(:,:,c)=subStruct.ITCsig.(chanFieldName);
            childPowBasesig(:,:,c)=subStruct.PowBaseSig.(chanFieldName);
            c=c+1;
        end
        
    end
    time = adultTimes(:,:,1); % times and freqs are same for children and adults
    freq = adultFreqs(:,:,1);

    % average channel data for adults and children
    avgERSPadult = mean(adultERSP,3);
    avgITCadult = mean(adultITC,3);
    avgPowBaseadult = mean(adultPowBase,3);
    avgERSPsigadult = mean(adultERSPsig,3);
    avgITCsigadult = mean(adultITCsig,3);
    avgPowBasesigadult = mean(adultPowBasesig,3);
    avgERSPchild = mean(childERSP,3);
    avgITCchild = mean(childITC,3);
    avgPowBasechild = mean(childPowBase,3);
    avgERSPsigchild = mean(childERSPsig,3);
    avgITCsigchild = mean(childITCsig,3);
    avgPowBasesigchild = mean(childPowBasesig,3);
    % saving average channel data into one file
    save(fullfile(adultDataPath,adultAvgName),'avgERSPadult','avgITCadult','avgPowBaseadult',...
        'avgERSPsigadult','avgITCsigadult','avgPowBasesigadult','time','freq')
    save(fullfile(childDataPath,childAvgName),"avgERSPchild","avgITCchild","avgPowBasechild",...
        "avgERSPsigchild","avgITCsigchild","avgPowBasesigchild","time","freq")

    % plotting channel average for children and adults
    % adults
    figA = figure;
    title1 = sprintf('Adult ERSP Channel %s',chanFieldName);
    title2 = sprintf('Adult ITC Channel %s',chanFieldName);
    subplot(2,1,1)
    tftopo(avgERSPadult,time,freq,'title',title1,'axcopy','off')
    subplot(2,1,2)
    tftopo(avgITCadult,time,freq,'title',title2,'axcopy','off')
    
    % children
    figB = figure;
    title3 = sprintf('Children ERSP Channel %s',chanFieldName);
    title4 = sprintf('Children ITC Channel %s',chanFieldName);
    subplot(2,1,1)
    tftopo(avgERSPchild,time,freq,'title',title3,'axcopy','off')
    subplot(2,1,2)
    tftopo(avgITCchild,time,freq,'title',title4,'axcopy','off')

    % saving figures
    adultFigName = sprintf('AdultAverages%s',chanFieldName);
    childFigName = sprintf('ChildAverages%s',chanFieldName);
    savefig(figA,fullfile(adultDataPath,adultFigName))
    savefig(figB,fullfile(childDataPath,childFigName))
    close all
    clear figA figB
end
% figure;
% 
% heatmap(avgERSPadult,times,freq,"Colormap",parula,"ColorbarVisible",'on','GridVisible','off')
% heatmap(time,freq,avgERSPadult')
% dfag=avgERSPadult';