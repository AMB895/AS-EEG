% 10/31/2024
% loading anti scores
clear
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
trialDataASpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/';
merged7t = readtable('merged_7t.csv'); % to get age and sex
idAgeDateSex = {merged7t.lunaid, merged7t.eeg_date, merged7t.eeg_age, merged7t.sex};
noEEGscan = find(isnan(cell2mat(idAgeDateSex(:,3))));
idAgeDateSex{1,1}(noEEGscan) = [];
idAgeDateSex{1,2}(noEEGscan) = [];
idAgeDateSex{1,3}(noEEGscan) = [];
idAgeDateSex{1,4}(noEEGscan) = [];
contents = dir(trialDataASpath);

%% Calculating Error Rate and average Latency
tableSize = [length(contents)-2 10];
varTypes = {'double','double','double','string','double','double','double','double','double','double'};
varNames = {'Luna ID','Scan Date','Age','Sex','Error Rate (no dropped)','Average Latency (correct trials)','# Correct','# Incorrect','# Dropped','Total Trials'};
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
    index_merged7t = find(idAgeDateSex{:,1} == subjID & idAgeDateSex{:,2} == scanDate);
    ageAtScan = idAgeDateSex{1,3}(index_merged7t);
    sex = string(idAgeDateSex{1,4}(index_merged7t));

    % Some participants are not in merge7t but have AS EEG score
    if isempty(ageAtScan)
        fprintf('Subject %d Date %d\n',subjID,scanDate)
         continue
    end
    % Loading AS score data
    subjectData = load(fileName);
    
    % Type of trial index
    asScore = subjectData.thisdata(:,6);
    corTrialNum = find(asScore == 1);
    incorTrialNum = find(asScore == 0);
    droppedTrialNum = find(isnan(asScore));

    % Calculating Average Latency for correct trials
    latency = subjectData.thisdata(:,5);
    latency(incorTrialNum) = [];
    % for d = 1:length(incorTrialNum) % Removing latency of incorrect trials
    %     latency(incorTrialNum(d)) = [];
    % end
    latency = latency(~isinf(latency)); % Removing latency of dropped trials
    averageLatency = mean(latency);

    % Error Rate Calculation
    numCorrect = sum(asScore == 1);
    numIncorrect = sum(asScore == 0);
    numDropped = length(find(isnan(asScore)));
    totalTrials = numDropped + numIncorrect + numCorrect;
    ER_noDropped = numIncorrect/(numIncorrect+numCorrect);

    ErrorLatencyTable(i-2,:) = {subjID,scanDate,ageAtScan,sex,ER_noDropped,averageLatency,numCorrect,numIncorrect,numDropped,totalTrials};
    trialType(i-2,:) = {subjID,scanDate,corTrialNum,incorTrialNum,droppedTrialNum};
end
idxMissingData = find(ErrorLatencyTable.("Luna ID") == 0);
ErrorLatencyTable(idxMissingData,:)=[];
trialTypeTable = array2table(trialType,'VariableNames',{'Luna ID','Scan Date','Correct Trials','Incorrect Trials','Dropped Trials'});


%% Save latency and error rate as .csv file
save('ErrorLatencyTable.mat','ErrorLatencyTable')
save('TypeOfTrial.mat','trialTypeTable')
writetable(ErrorLatencyTable,'7t_eegAS_ErrorLatency.csv')