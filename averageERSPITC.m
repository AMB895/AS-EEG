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
% Getting count of total adults and children

for i=1:size(allSubFolders)
    nameDate = allSubFolders(i).name;
    splitNameDate = split(nameDate,'_');
    name(i) = str2double(cell2mat(splitNameDate(1)));
    date(i) = str2double(cell2mat(splitNameDate(2)));
    ageIdx = find(ErrorLatencyTable.("Luna ID")==name(i)&ErrorLatencyTable.("Scan Date")==date(i));
    if isempty(ageIdx)
        continue;
    end
    AGE(i) = ErrorLatencyTable.Age(ageIdx);
end
nameDateAge=[name',date',AGE'];
numAdults = 0;
numChildren = 0;
for j=1:length(nameDateAge)
    if nameDateAge(j,3) >= 18
        numAdults = numAdults+1;
    elseif nameDateAge(j,3) < 18
        numChildren = numChildren +1;
    end
end
%%
% Splitting up children and adults
for chan=1:64 
    chanFieldName= sprintf('Ch%d',chan);
    % names of files being saved
    adultAvgName = sprintf('AdultAverages%s.mat',chanFieldName);
    childAvgName = sprintf('ChildAverages%s.mat',chanFieldName);

    % check to see if channel data has already been created
    if exist(fullfile(childDataPath,childAvgName),'file') && exist(fullfile(adultDataPath,adultAvgName),'file')
        fprintf('skipping; already created %s and %s\m',adultAvgName,childAvgName)
        continue;
    end

    fprintf('\n-------Processing %s-------\n',chanFieldName)
    % defining intial indicies for adult and child matricies
    a=1;
    c=1;
    for currentSub = 1:size(allSubFolders)
        % Getting count of adults and children
         
        subjectFolderName = sprintf('%s',allSubFolders(currentSub).name);
        structureName = sprintf('%s.mat',allSubFolders(currentSub).name);
        subFieldName = sprintf('sub%s',subjectFolderName);
        % loaded variable is called subStruct
        load(sprintf('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/TimeFreqPlots/%s/%s',subjectFolderName,structureName));
        age = subStruct.Age;
        if isempty(age)
            fprintf('No age %s\n',subjectFolderName)
        end
        if age >= 18
            % do these have to be cell arrays?
            adultERSP(:,:,a)=subStruct.ERSP.(chanFieldName); % 3rd dimension is subject
            adultITC(:,:,a)=subStruct.ITC.(chanFieldName);
            adultTimes(:,:,a)=subStruct.Times.(chanFieldName);
            adultFreqs(:,:,a)=subStruct.Freqs.(chanFieldName);
            adultPowBase(:,:,a)=subStruct.PowBase.(chanFieldName);
            adultERSPboot(:,:,a) = subStruct.ERSPbootsig.(chanFieldName);
            adultITCboot(:,:,a) = subStruct.ITCbootsig.(chanFieldName);
            a=a+1;
        elseif age < 18
            childERSP(:,:,c)=subStruct.ERSP.(chanFieldName);
            childITC(:,:,c)=subStruct.ITC.(chanFieldName);
            childTimes(:,:,c)=subStruct.Times.(chanFieldName);
            childFreqs(:,:,c)=subStruct.Freqs.(chanFieldName);
            childPowBase(:,:,c)=subStruct.PowBase.(chanFieldName);
            childERSPboot(:,:,c) = subStruct.ERSPbootsig.(chanFieldName);
            childITCboot(:,:,c) = subStruct.ITCbootsig.(chanFieldName);
            c=c+1;
        end
        
    end
    fprintf('Processed %d Adults\n',(a-1))
    fprintf('Processed %d Children\n',(c-1))
    time = adultTimes(:,:,1); % times and freqs are same for children and adults
    freq = adultFreqs(:,:,1);

    % average channel data for adults
    avgERSPadult = mean(adultERSP,3);
    avgITCadult = mean(adultITC,3);
    avgPowBaseadult = mean(adultPowBase,3);
    avgERSPbootadult = mean(adultERSPboot,3);
    avgITCbootadult = mean(adultITCboot,3);

    % average channel data for children
    avgERSPchild = mean(childERSP,3);
    avgITCchild = mean(childITC,3);
    avgPowBasechild = mean(childPowBase,3);
    avgERSPbootchild = mean(childERSPboot,3);
    avgITCbootchild = mean(childITCboot,3);

    % saving average channel data into one file
    save(fullfile(adultDataPath,adultAvgName),'avgERSPadult','avgITCadult','avgPowBaseadult',...
       'avgERSPbootadult','avgITCbootadult','time','freq')
    save(fullfile(childDataPath,childAvgName),"avgERSPchild","avgITCchild","avgPowBasechild",...
        "avgERSPbootchild","avgITCbootchild","time","freq")

    % plotting channel average for children and adults
    % adults
    figA(1) = figure;
    % no significance mask
    subplot(2,1,1)
    title1 = sprintf('Adult ERSP %s',chanFieldName);
    tftopo(avgERSPadult,time,freq,'title',title1,'axcopy','off')
    colorbar
    subplot(2,1,2)
    title2 = sprintf('Adult ITC %s',chanFieldName);
    tftopo(avgITCadult,time,freq,'title',title2,'axcopy','off')
    colorbar
    % significance mask- only plotting ERSP mask for now, ITC isn't working because it is complex
    figA(2) = figure;
    % subplot(2,1,1)
    title3 = sprintf('Adult Significant ERSP %s',chanFieldName);
    tftopo(avgERSPadult,time,freq,'signifs',avgERSPbootadult,'title',title3,'axcopy','off')
    colorbar
    % subplot(2,1,2)
    % title4 = sprintf('Adult Significant ITC %s',chanFieldName);
    % tftopo(avgITCadult,time,freq,'signifs',avgITCbootadult,'title',title4,'axcopy','off')
    % colorbar

    % children
    figC(1) = figure;
    % no significance mask
    subplot(2,1,1)
    title5 = sprintf('Children ERSP %s',chanFieldName);
    tftopo(avgERSPchild,time,freq,'title',title5,'axcopy','off')
    colorbar
    subplot(2,1,2)
    title6 = sprintf('Children ITC %s',chanFieldName);
    tftopo(avgITCchild,time,freq,'title',title6,'axcopy','off')
    colorbar
    % significance mask - only plotting ERSP mask for now, ITC isn't working because it is complex
    figC(2) = figure;
    % subplot(2,1,1)
    title7 = sprintf('Children Significant ERSP %s',chanFieldName);
    tftopo(avgERSPchild,time,freq,'signifs',avgERSPbootchild,'title',title7,'axcopy','off')
    colorbar
    % subplot(2,1,2)
    % title8 = sprintf('Children Significant ITC %s',chanFieldName);
    % tftopo(avgITCchild,time,freq,'signifs',avgITCbootchild,'title',title8,'axcopy','off')
    % colorbar

    % saving figures
    adultFigName = sprintf('AdultAverages%s',chanFieldName);
    childFigName = sprintf('ChildAverages%s',chanFieldName);
    savefig(figA,fullfile(adultDataPath,adultFigName))
    savefig(figC,fullfile(childDataPath,childFigName))
    close all
    clear figA figC
end

%% only looking at plots between 3-30Hz and 0-500 ms
% messing around
close all
timeIdx0 = find(time==0);
timeIdx500 = find(time==500);
newTime = time(timeIdx0:timeIdx500);
freqIdx3 = find(freq==3);
freqIdx30 = find(freq>30,1,"first");
newFreq = freq(freqIdx3:freqIdx30);
newERSP = avgERSPadult(freqIdx3:freqIdx30,timeIdx0:timeIdx500);
figure;
surf(newTime,newFreq,newERSP)
view(2)
title('adults')
colorbar
figure;
surf(newTime,newFreq,avgERSPchild(freqIdx3:freqIdx30,timeIdx0:timeIdx500))
view(2)
title('children')
colorbar