% Abby Beatty
% 10/31/2024
% Updated 1/30/25 bc now have error corrected trials
% Updataed 2/18/25 
% Updated 3/20/25 to calculate coefficient of latency variance and percent error corrected
% TODO: add scoring script to github for this project

clear; close all;
% adding paths
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/') % directory to EOG calibration
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/') % directory to EOG trial level data
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/') % for merge 7t
trialDataASpath = '/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/1*_20*.mat'; % file paths for trial level data

% loading AS trial data from score_anti() 
antiScoreTable = readtable('eeg_anti_20250802.csv');
merged7t = readtable('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/merged_7t.csv'); % to get age and sex

% only care about lunaid, eeg scan date, age at eeg scan and sex from merged7t
merged7t_eeg = [merged7t(:,'lunaid'), merged7t(:,'eeg_date'), merged7t(:,'eeg_age'), merged7t(:,'sex')];

% Total number of AS eegs
contents = dir(trialDataASpath);

totASeegs = size(contents,1);

%% Calculating Percent error corrected trials, average correct latency and coefficent of latency variation
% (% error corrected) = (# errcor) / (#errcor + #cor)
% Coeff of variation = (subject's latency std)/(subject's avg correct lat)

% Setting Error Latency table up
varTypes = {'double','double','double','string','double','double','double','double','double','double','double','double','double','double','double','double'};
varNames = {'LunaID','ScanDate','Age','Sex','PercentErrCor','Accuracy','AvgLatCor','CVLatCor','AvgLatErrCor','CVLatErrCor','NumCor','NumErrCor','NumIncor','NumDrop','Total','Viable'};
tableSize = [length(contents) length(varNames)];
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
    
    % Total number of correct, error corrected, incorrect and droppped trials
    numCor = length(cor_idx);
    numErrCor = length(errcor_idx);
    numIncor = length(incor_idx);
    numDrop = length(drop_idx);
    totalTrials = numCor + numErrCor + numIncor + numDrop;

    % Latency for every correct and error corrected trial
    latency = subTable.Latency;
    corLats = latency(cor_idx);
    errcorLats = latency(errcor_idx);
        
    % average latency for correct and error corrected trials
    avgCorLat = mean(corLats);   
    avgErrCorLat = mean(errcorLats);
    
    % SD of latency for correct and error corrected trials
    SDcorLat = std(corLats);
    SDerrcorLat = std(errcorLats);
    
    % Coefficient of Latency Variation for correct and error corrected trials
    cvLatCor = SDcorLat/avgCorLat;
    cvLatErrCor = SDerrcorLat/avgErrCorLat;
    
    % Percent Error Corrected Calculation
    percentErrCor = numErrCor / (numErrCor + numCor);
    
    % Accuracy Calculation
    accuracy = numCor / (numErrCor + numCor);
    
    % Determining if subject is viable to use for stats
    numOnTask = numCor + numErrCor;
    if numOnTask >= totalTrials*0.5
        viable = 1;
    elseif numOnTask < totalTrials*0.5
        viable = 0;
    end
        
    ErrorLatencyTable(currentSub,:) = {lunaID,scanDate,age,sex,percentErrCor,accuracy,avgCorLat,cvLatCor,avgErrCorLat,cvLatErrCor,numCor,numErrCor,numIncor,numDrop,totalTrials,viable};
end

idxMissingData = find(ErrorLatencyTable.LunaID == 0);
ErrorLatencyTable(idxMissingData,:)=[];

%% Save latency and error rate as .csv file
% get current date
currentdate = datetime();
currentdate.Format = 'yyyyMMdd';
savename_csv = sprintf('eegAS_ErrorLatency_%s.csv',currentdate);
savename_mat = sprintf('eegAS_ErrorLatency_%s.mat',currentdate);
save(sprintf('/Volumes/Hera/Abby/AS_EEG/Results/%s',savename_mat),'ErrorLatencyTable')
writetable(ErrorLatencyTable,sprintf('/Volumes/Hera/Abby/AS_EEG/%s',savename_csv))
