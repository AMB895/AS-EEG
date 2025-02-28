%% Stats on ERSP data per condition (correct, incorrect, error corrected)
% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

% Uncomment what type of trial you want to run
trialtype = 1; % correct trials
% trialtype = 0; % incorrect trials
% trialtype = 2; % error corrected trials

% Setting up directories and paths
if trialtype ==1
    % all ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials';
    allerspdataname = 'all_ERSP_DATA_cor.mat';
    % Files being made
    ttestdataname = 'ttest_outputs_cor.mat';
    ageregressdataname = 'linregress_age_outputs_cor.mat';
    invageregressdataname = 'linregress_invage_outputs_cor.mat';
elseif trialtype ==0
    % all ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/IncorrectTrials';
    allerspdataname = 'all_ERSP_DATA_incor.mat';
    % Files being made
    ttestdataname = 'ttest_outputs_incor.mat';
    ageregressdataname = 'linregress_age_outputs_incor.mat';
    invageregressdataname = 'linregress_invage_outputs_incor.mat';
elseif trialtype ==2
    % all ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials';
    allerspdataname = 'all_ERSP_DATA_errcor.mat';
    % Files being made
    ttestdataname = 'ttest_outputs_errcor.mat';
    ageregressdataname = 'linregress_age_outputs_errcor.mat';
    invageregressdataname = 'linregress_invage_outputs_errcor.mat';
end
%% One sample T-test for areas of activation (group level)
% One-sample t-test for every channel/time/freq
% Group activity significantly greater than baseline

% load ERSP data and ID info
% erspdata: [subject,chan, freq, time]
load(fullfile(allerspdatapath,allerspdataname))
load(fullfile(allerspdatapath,'ID_info.mat'))
erspdata = allerspdata;

% defining number of channels, frequencies, and times to loop through
numChans = size(erspdata,2);
numFreqs = size(erspdata,3);
numTimes = size(erspdata,4);

% skipping if group activation t-test outputs exist
if exist(fullfile(allerspdatapath,ttestdataname),'file')
    fprintf('skipping; Already computed one-sample t-test\n')
else
    % setting up progress bar
    f = waitbar(0,'Chan','Name','One-sample T-test Channel Progress');
    
    for i = 1:numChans
        % updating progress bar
        waitbar(i/numChans,f,sprintf('Channel %d/%d',i,numChans))
        for j = 1:numFreqs
            for k = 1:numTimes
                % t-test for every (time,frequency) point for each channel
                data = squeeze(erspdata(:,i,j,k));
                % saving t-statistic and p-value
                [~,pval_ttest(i,j,k),~,ttest_stats] = ttest(data);
                tval(i,j,k) = ttest_stats.tstat;
            end
        end
    end
    save(fullfile(allerspdatapath,ttestdataname),'pval_ttest','tval')
end

%% Age effects for channel/frequency/time- Linear Regression
% erspdata: [subject,chan,freq,time]
load(fullfile(allerspdatapath,allerspdataname))
load(fullfile(allerspdatapath,'ID_info.mat'))
erspdata = allerspdata;

% skipping if age effects linear regression outputs exist
if exist(fullfile(allerspdatapath,ageregressdataname),'file')
    fprintf('skipping; Already computed linear regression for age-effects\n')
else
    
    % defining number of channels, frequencies, and times to loop through
    numChans = size(erspdata,2);
    numFreqs = size(erspdata,3);
    numTimes = size(erspdata,4);
    
    % Making age vector to input into regress()
    % regress(y,x) inputs:
    %   y = response vector
    %   x = predictor vector w/ column of ones to calculate coefficients
    age = [ones(size(erspdata,1),1),IDmatrix(:,3)];
    
    % setting up progress bar
    f = waitbar(0,'Chan','Name','Linear Regression Progress');
    
    for i = 1:numChans
        % updating progress bar
        waitbar(i/numChans,f,sprintf('Channel %d/%d',i,numChans))
        for j = 1:numFreqs
            for k = 1:numTimes
                % Linear regression for every (time,freq) point for every channel
                data = squeeze(erspdata(:,i,j,k));
                [b,~,~,~,regress_stats] = regress(data,age);
                % regress_stats : [R-square, F stat, p value, Error variance]
                % saving b(2) (age coefficient, b(1) = coefficient of intercept)
                b_age(i,j,k) = b(2);
                % saving f-statistic and corresponding p-value
                fval(i,j,k) = regress_stats(1,2);
                pval(i,j,k) = regress_stats(1,3);
            end
        end
    end
    save(fullfile(allerspdatapath,ageregressdataname),'b_age','fval','pval')
end

%% Inverse Age effects for channel/frequency/time- Linear Regression
% erspdata: [subject,chan,freq,time]
load(fullfile(allerspdatapath,allerspdataname))
load(fullfile(allerspdatapath,'ID_info.mat'))
erspdata = allerspdata;

% skipping if inverse age effects linear regression outputs exist
if exist(fullfile(allerspdatapath,invageregressdataname),'file')
    fprintf('skipping; Already computed linear regression for inverse age-effects\n')
else

    % defining number of channels, frequencies, and times to loop through
    numChans = size(erspdata,2);
    numFreqs = size(erspdata,3);
    numTimes = size(erspdata,4);

    % Making age vector to input into regress()
    % regress(y,x) inputs:
    %   y = response vector
    %   x = predictor vector w/ column of ones to calculate coefficients
    inverse_age = [ones(size(erspdata,1),1),(1./IDmatrix(:,3))];
    
    % setting up progress bar
    f = waitbar(0,'Chan','Name','Inverse Age Linear Regression Progress');
    
    for i = 1:numChans
        % updating progress bar
        waitbar(i/numChans,f,sprintf('Channel %d/%d',i,numChans))
        for j = 1:numFreqs
            for k = 1:numTimes
                % Linear regression for every (time,freq) point for every channel
                data = squeeze(erspdata(:,i,j,k));
                [b_inv,~,~,~,regress_stats_inv] = regress(data,inverse_age);
                % regress_stats : [R-square, F stat, p value, Error variance]
                % saving b(2) (age coefficient, b(1) = coefficient of intercept)
                b_age_inv(i,j,k) = b_inv(2);
                % saving f-statistic and corresponding p-value
                fval_inv(i,j,k) = regress_stats_inv(1,2);
                pval_inv(i,j,k) = regress_stats_inv(1,3);
            end
        end
    end
    save(fullfile(allerspdatapath,invageregressdataname),'b_age_inv','fval_inv','pval_inv')
end