% 10/31/2024
% loading anti scores
clear
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
trialDataASpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/';
merged7t = readtable('merged_7t.csv'); % to get ages
idAgeDate = [merged7t.lunaid, merged7t.eeg_date, merged7t.eeg_age];
noEEGscan = find(isnan(idAgeDate(:,3)));
contents = dir(trialDataASpath);

%% Calculating Error Rate and average Latency
tableSize = [length(contents)-2 9];
varTypes = {'string','string','double','double','double','double','double','double','double'};
varNames = {'Luna ID','Scan Date','Age','Error Rate (no dropped)','Average Latency (correct trials)','# Correct','# Incorrect','# Dropped','Total Trials'};
ErrorLatencyTable = table('Size',tableSize,'VariableTypes',varTypes,'VariableNames',varNames);
trialType = cell([length(contents)-2 5]);

for i=1:length(contents)
    if contents(i).bytes == 0 % first two fields in contents are not trial data
        continue
    end
    fileName = contents(i).name;
    nameDate = strsplit(fileName,{'_','.'});
    subjID = str2double(nameDate{1});
    scanDate = str2double(nameDate{2});
    index_merged7t = find(idAgeDate(:,1) == subjID & idAgeDate(:,2) == scanDate);
    ageAtScan = idAgeDate(index_merged7t,3);
    if isempty(ageAtScan)
        ageAtScan = NaN;
    end
    subjectData = load(fileName);
    asScore = subjectData.thisdata(:,6);
    latencyWithInfs = subjectData.thisdata(:,5);
    infsLocation = isinf(latencyWithInfs);
    latency = latencyWithInfs(~infsLocation);
    averageLatency = mean(latency);
    numCorrect = sum(asScore == 1);
    corTrialNum = find(asScore == 1);
    numIncorrect = sum(asScore == 0);
    incorTrialNum = find(asScore == 0);
    numDropped = length(find(isnan(asScore)));
    droppedTrialNum = find(isnan(asScore));
    totalTrials = numDropped + numIncorrect + numCorrect;
    ER_noDropped = numIncorrect/(numIncorrect+numCorrect);
    % ER_withDropped = (numIncorrect+numDropped)/(totalTrials);
    ErrorLatencyTable(i-2,:) = {subjID,scanDate,ageAtScan,ER_noDropped,averageLatency,numCorrect,numIncorrect,numDropped,totalTrials};
    trialType(i-2,:) = {subjID,scanDate,corTrialNum,incorTrialNum,droppedTrialNum};
end
trialTypeTable = array2table(trialType,'VariableNames',{'Luna ID','Scan Date','Correct Trials','Incorrect Trials','Dropped Trials'});


%% Save latency and error rate as .csv file
save('ErrorLatencyTable.mat','ErrorLatencyTable')
save('TypeOfTrial.mat','trialTypeTable')
writetable(ErrorLatencyTable,'7t_eegAS_ErrorLatency.csv')