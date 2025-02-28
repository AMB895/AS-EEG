% 10/31/2024
% Updated 1/30/25 bc now have error corrected trials
% Updataed 2/18/25 

% adding paths
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/') % for merge 7t
trialDataASpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/1*_20*.mat';

% loading AS trial data from score_anti()
antiScoreTable = readtable('eeg_anti_20250218.csv');
merged7t=readtable('merged_7t.csv'); % to get age and sex

% only care about lunaid, eeg scan date, age at eeg scan and sex from merged7t
merged7t_eeg = [merged7t(:,'lunaid'), merged7t(:,'eeg_date'), merged7t(:,'eeg_age'), merged7t(:,'sex')];

% finding subjects who did not do not have an eeg scan and removing them from merge_7t_eeg
noeeg_idx = find(isnan(merged7t_eeg.eeg_date));
merged7t_eeg(noeeg_idx,:) = [];

contents = dir(trialDataASpath);


%% Calculating Error Rate and average Latency
% Bea error rate = # err/(# cor + # err cor)
% error rate = # err / (# cor + # err + # err cor)

% Setting Error Latency table up
varTypes = {'double','double','double','string','double','double','double','double','double','double','double','double','double'};
varNames = {'LunaID','ScanDate','Age','Sex','ErrorRateBea','ErrorRate','AverageLatencyCorrect','Correct','ErrorCorrected','Incorrect','Dropped','Total','Viable'};
tableSize = [length(contents) length(varTypes)];
ErrorLatencyTable = table('Size',tableSize,'VariableTypes',varTypes,'VariableNames',varNames);

missingMerged7t = [];
missingASscore = [];
for currentSub=1:length(contents)
    
    fileName = contents(currentSub).name;
    nameDate = strsplit(fileName,{'_','.'});
    lunaID = str2double(nameDate{1});
    scanDate = str2double(nameDate{2});
    idx_merged7t = find(merged7t_eeg.lunaid == lunaID & merged7t_eeg.eeg_date == scanDate);
    age = merged7t_eeg.eeg_age(idx_merged7t);
    sex = string(merged7t_eeg.sex(idx_merged7t));

    % Some participants are not in merge7t (have no age) but have AS EEG score
    if isempty(age)
        fprintf('\nSubject %d Date %d is not in merged7t\n',lunaID,scanDate)
        missingMerged7t(end+1,:) = [lunaID,scanDate];
        continue
    end
    
    % Loading AS score data
    scoreTable_idx = find(antiScoreTable.LunaID == lunaID & antiScoreTable.ScanDate == scanDate);
    subTable = antiScoreTable(scoreTable_idx,:);
        
    % some participants are in merged7t but do not have AS eeg score
    if isempty(scoreTable_idx)
       fprintf('\nSubject %d Date %d does not have AS eeg scored\n',lunaID,scanDate)
       missingASscore(end+1,:) = [lunaID,scanDate];
       continue
    end
    
    % Type of trial index
    cor_idx = find(subTable.Correct == 1);
    incor_idx = find(subTable.Correct == 0);
    errcor_idx = find(subTable.Correct == 2);
    drop_idx = find(subTable.Correct == -1);

    % Calculating Average Latency for correct trials
    latency = subTable.Latency;
    allRemoveLat_idx = [incor_idx; errcor_idx; drop_idx];
    latency(allRemoveLat_idx) = [];
        
    % calculating average latency for correct trials
    avgLat = mean(latency);   

    % Error Rate Calculation
    numCor = length(cor_idx);
    numIncor = length(incor_idx);
    numErrCor = length(errcor_idx);
    numDrop = length(drop_idx);
    ER_bea = numIncor/(numErrCor + numCor);
    ER = numIncor / (numCor + numIncor + numErrCor);
    totalTrials = numCor + numIncor + numErrCor + numDrop;
    
    % Determining if subject is viable to use for stats
    numViable = numCor+numIncor+numErrCor;
    numOnTask = numCor + numErrCor;
    if numOnTask >= totalTrials*0.5
        viable = 1;
    elseif numOnTask < totalTrials*0.5
        viable = 0;
    end
        
    ErrorLatencyTable(end+1,:) = {lunaID,scanDate,age,sex,ER_bea,ER,avgLat,numCor,numErrCor,numIncor,numDrop,totalTrials,viable};
end
idxMissingData = find(ErrorLatencyTable.LunaID == 0);
ErrorLatencyTable(idxMissingData,:)=[];

%% Save latency and error rate as .csv file
save('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable_20250218.mat','ErrorLatencyTable')
writetable(ErrorLatencyTable,'/Volumes/Hera/Abby/AS_EEG/7t_eegAS_ErrorLatency_20250218.csv')

%% Stats on latency and error rate
% Viable subjects table
viableSubsIdx = find(ErrorLatencyTable.Viable ==1);
viableTable = ErrorLatencyTable(viableSubsIdx,:);

% Latency- Linear Age
age = [ones(size(viableTable,1),1),viableTable.Age];
lat = viableTable.AverageLatencyCorrect;
[b_lat,~,~,~,stats_lat] = regress(lat,age);
% stats: [R squared, F stat, p-value,error variance]

% Latency- Inverse Age
invage = [ones(size(viableTable,1),1),1./viableTable.Age];
[b_lat_inv,~,~,~,stats_lat_inv] = regress(lat,invage);

% Error Rate- Linear Age
er = viableTable.ErrorRate;
[b_er,~,~,~,stats_er] = regress(er,age);

% Error Rate- Inverse Age
[b_er_inv,~,~,~,stats_er_inv] = regress(er,invage);
