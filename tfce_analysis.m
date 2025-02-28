%% TFCE on t/F values from t-test and linear regression
clear; close all;
%% Setting up
% increase number of cores for parallel processing
% limo defaults to 35, but only see 22
setenv('NUMBER_OF_PROCESSORS','50')

% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

% Setting up directories and paths
% TFCE outputs data path
tfcepath = '/Volumes/Hera/Abby/AS_EEG/TFCE/';
cor_origtfcename = 'TFCE_scores_correct.mat';
cor_incor_origtfcename = 'TFCE_scores_cor_incor.mat';
cor_errcor_origtfcename = 'TFCE_scores_cor_errcor.mat';

% load t-test t-values for all trial types
cor_ttest = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/CorrectTrials/ttest_outputs_cor.mat');
cor_tval = cor_ttest.tval;

% load two-sample t-test t-values
cor_incor_ttest = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/two_sample_ttest_cor_incor.mat');
cor_incor_tval = cor_incor_ttest.t;
cor_errcor_ttest = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/two_sample_ttest_cor_errcor.mat');
cor_errcor_tval = cor_errcor_ttest.t;

% load linear regression with age F-values for al trial types
cor_age = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/CorrectTrials/linregress_age_outputs_cor.mat');
cor_fval_age = cor_age.fval;

% load linear regression with inverse age F-values for al trial types
cor_invage = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/CorrectTrials/linregress_invage_outputs_cor.mat');
cor_fval_invage = cor_invage.fval_inv;

% load channel logical matrix
load('/Volumes/Hera/Abby/AS_EEG/STUDY/logicalchanneighbmatrix.mat')

%% Run TFCE on original t/F stats for correct trials
% Inputs to limo_tfce():
%   type: 1,2, or 3 for 1D, 2D, or 3D data
%       for 3D data (channels x freqs x time)
%   data: map of t/F values
%   channeighbstructmat: neighborhood matrix for clustering

% skip if completed for correct trials
if exist(fullfile(tfcepath,cor_origtfcename),'file') 
    fprintf('skipping; already computed TFCE scores on original data from correct trials\n')
else
    % Group activation: t-test
    % could run on lots of levels of analyses
    tfce_scores_group_cor = limo_tfce(3,cor_tval,channeighbmatrix);
    % Age effects: Linear Regression
    tfce_scores_age_cor = limo_tfce(3,cor_fval_age,channeighbmatrix);
    tfce_scores_ageinv_cor = limo_tfce(3,cor_fval_invage,channeighbmatrix);

    % save original TFCE scores
    save(fullfile(tfcepath,cor_origtfcename),'tfce_scores_group_cor','tfce_scores_age_cor','tfce_scores_ageinv_cor')
end

%% Run TFCE on t-values from two-sample t-test between conditions
% Correct and Incorrect
if exist(fullfile(tfcepath,cor_incor_origtfcename),'file')
    fprintf('Completed TFCE on two-sample t-test between correct and incorrect trials\n')
else
    tfce_scores_group_cor_incor = limo_tfce(3,cor_incor_tval,channeighbmatrix);
    
    % save original TFCE socres
    save(fullfile(tfcepath,cor_incor_origtfcename),'tfce_scores_group_cor_incor')
end

% Correct and Error Corrected
if exist(fullfile(tfcepath,cor_errcor_origtfcename),'file')
    fprintf('Completed TFCE on two-sample t-test between correct and error corrected trials\n')
else
    tfce_scores_group_cor_errcor = limo_tfce(3,cor_errcor_tval,channeighbmatrix);
    
    % save original TFCE socres
    save(fullfile(tfcepath,cor_errcor_origtfcename),'tfce_scores_group_cor_errcor')
end