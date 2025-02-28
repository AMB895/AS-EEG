function [EEG] = remarkAS(taskdirectory,dryrun)
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal') % need to run eeg_data.m

path_file = hera('Raw/EEG/7TBrainMech');
%path_file = '/Volumes/Hera/Raw/EEG/7TBrainMech';
outputpath = [taskdirectory '/remarked/'];

% directory of EEG data
[path ,folder] = fileparts(path_file);
d =[path,'/',folder,'/'];

namesOri = dir([d,'*/*.bdf']);

antiIDX = find (cellfun (@any,regexpi ( {namesOri.name}.', 'anti')));

for idx = antiIDX'
    
    currentName = namesOri(idx).name(1:end-4);
    currentNameSplit = split(currentName,'_');
    subjectID = [currentNameSplit{1} '_' currentNameSplit{2}]; % need for eeg_data.m
    d = [namesOri(idx).folder '/'];
    
    % % skipping 11834_20210628 b/c anti trial is 1545 secs (most are ~120 secs)
    if currentName == "11834_20210628_anti"
        continue
    end

    % skip if already been remarked
    finalfile=fullfile(outputpath, [currentName '_Rem.set']);
     if exist(finalfile,'file')
         fprintf('already have %s\n', finalfile)
         continue
     end
     if dryrun
         fprintf('want to run %s; set dryrun=0 to actually run\n', finalfile)
         continue
     end
    fprintf('making %s\n',finalfile);
    
    % loading EEG data
    f = [d currentName '.bdf'];
    EEG = pop_biosig(f);
    EEG.setname=[currentName 'Rem']; %name the EEGLAB setclc (this is not the set file itself)
    eegData = eeg_data('#anti',{'Status'},'subjs',{subjectID}); % loading fixed status channel
    triggers = unique(eegData.Status); 
    
    % Trigger Values:
    % 254 = ITI
    % 101-105 = Fixation
    % 151-155 = Anti Target
    
    itiLoc = find(eegData.Status == 254); % Find location of ITI in raw status channel
    itiDiff = diff(itiLoc); % Don't want locations from the same trial
    itiDiffIndex = [1 (find(itiDiff > 1)+1)]; % Need to manually force location 1
    itiOnset = itiLoc(itiDiffIndex); % Only want iti at the start of each trial

    fixationLoc = find(eegData.Status > 100 & eegData.Status < 106); % Find location of fixation in raw status channel
    fixationDiff = diff(fixationLoc); % Don't want locations from the same trial
    fixationDiffIndex = [1 (find(fixationDiff > 1) +1)]; % Need to manually force location 1
    fixationOnset = fixationLoc(fixationDiffIndex); % Only want fixation at the start of each trial

    targetLoc = find(eegData.Status > 150 & eegData.Status < 156); % Find location of target in raw status channel
    targetDiff = diff(targetLoc); % Don't want locations from the same trial
    targetDiffIndex = [1 (find(targetDiff >1)+1)]; % Need to manually force location 1
    targetOnset = targetLoc(targetDiffIndex); % Only want target at the start of each trial

% EEG.event.latency gives index of event in raw EEG signal
% Finding when the latency matches the indicies found above, changing event
% type when latency and index match
    eventOnsetIndex = cell2mat({EEG.event(:).latency}); 
    for i=1:length(itiOnset)
        itiLatencyIndex(i) = find(eventOnsetIndex == itiOnset(i));
        EEG = pop_editeventvals(EEG,'changefield',{itiLatencyIndex(i),'type',1});
    end

    for j=1:length(fixationOnset)
        fixationLatencyIndex(j) = find( eventOnsetIndex== fixationOnset(j));
        EEG = pop_editeventvals(EEG,'changefield',{fixationLatencyIndex(j),'type',2});
    end

    for k = 1:length(targetOnset)
        targetLatencyIndex(k) = find( eventOnsetIndex == targetOnset(k));
        EEG = pop_editeventvals(EEG,'changefield',{targetLatencyIndex(k),'type',3});
    end
    
    if ischar(EEG.event(1).type)
        for n = 1:length({EEG.event(:).type})
            EEG.event(n).type = str2num(EEG.event(n).type);
        end
    end
    eventTypesAfter = cell2mat({EEG.event(:).type}); % getting new event types after remarking
    
    % Some events do not correspond to ITI, fixation or target -> delete event
    for b = 1:length(eventTypesAfter)
        if eventTypesAfter(b) > 3
            EEG = pop_editeventvals(EEG,'delete',b);
        end
    end

    EEG = pop_saveset( EEG, 'filename',[currentName '_Rem.set'],'filepath',outputpath);
    
% checking to make sure remarking was done correctly (should have 40 of each mark)
% may get 41 ITI marks
    countItiMarks = sum(cell2mat({EEG.event(:).type})==1);
    countFixationMarks = sum(cell2mat({EEG.event(:).type})==2);
    countTargetMarks = sum(cell2mat({EEG.event(:).type})==3);
    fprintf('# ITI Marks: %d\n# Fixation Marks: %d\n# Target Marks: %d\n',countItiMarks,countFixationMarks,countTargetMarks)
    
    clear i j k
end