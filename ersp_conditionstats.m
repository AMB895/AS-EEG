%% Two-sample T-test for conditions on ERSP data
% Correct vs. Incorrect
% Correct vs. Error Corrected
% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat') % to get ages
% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

% Load data
cor_ersp = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/all_ERSP_DATA_cor.mat');
cor_id = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/ID_info.mat');
cor_id = cor_id.IDmatrix;
times = cor_ersp.times;
freqs = cor_ersp.freqs;
cor_ersp = cor_ersp.allerspdata;
errcor_ersp = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/all_ERSP_DATA_errcor.mat');
errcor_id = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/ID_info.mat');
errcor_ersp = errcor_ersp.allerspdata;
errcor_id = errcor_id.IDmatrix;
incor_ersp = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/IncorrectTrials/all_ERSP_DATA_incor.mat');
incor_id = load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/IncorrectTrials/ID_info.mat');
incor_ersp  = incor_ersp.allerspdata;
incor_id = incor_id.IDmatrix;

%% Correct vs. Incorrect Trials

% check to see if first 'subject' in allerspdata is empty
% for some reason in my ersp_analysis.m it added an empty ersp data frame
% to allerspdata
% correct trials
if (sum(squeeze(cor_ersp(:,:,:,1) ~=0),'all') == 0) && (length(cor_id) ~= size(cor_ersp,4))
    fprintf('First subject in correct ersp data is empty\nDeleting cor_ersp(:,:,:,1)\n')
    cor_ersp = cor_ersp(:,:,:,2:end);
    erspdata1 = cor_ersp;
else
    fprintf('No empty subject in cor_ersp\n')
    erspdata1 = cor_ersp;
end

% incorrect trials
if (sum(squeeze(incor_ersp(:,:,:,1) ~=0),'all') == 0) && (length(incor_id) ~= size(incor_ersp,4))
    fprintf('First subject in incorrect ersp data is empty\nDeleting incor_ersp(:,:,:,1)\n')
    incor_ersp = incor_ersp(:,:,:,2:end);
    erspdata2 = incor_ersp;
else
    fprintf('No empty subject in incor_ersp\n')
    erspdata2 = incor_ersp;
end

numChans = size(cor_ersp,1);
numFreqs = size(cor_ersp,2);
numTimes = size(cor_ersp,3);
if exist('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/two_sample_ttest_cor_incor.mat','file')
    fprintf('already computed two-sample t-test for correct and incorrect trials\n')
else
    % setting up progress bar
    f = waitbar(0,'Chan','Name','Two-sample T-test Channel Progress');
    for i=1:numChans
        % updating waitbar
        waitbar(i/numChans,f,sprintf('Channel %d/%d',i,numChans))
        for j = 1:numFreqs
            for k = 1:numTimes
                data1 = squeeze(erspdata1(i,j,k,:));
                data2 = squeeze(erspdata2(i,j,k,:));
                [~,p(i,j,k),~,stats] = ttest2(data1,data2);
                t(i,j,k) = stats.tstat;
            end
        end    
    end
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/two_sample_ttest_cor_incor.mat','p','t')
end

%% Correct vs. Error corrected Trials
clear p t

% check to see if first 'subject' in allerspdata is empty
% for some reason in my ersp_analysis.m it added an empty ersp data frame
% to allerspdata
% correct trials
if (sum(squeeze(cor_ersp(:,:,:,1) ~=0),'all') == 0) && (length(cor_id) ~= size(cor_ersp,4))
    fprintf('First subject in correct ersp data is empty\nDeleting cor_ersp(:,:,:,1)\n')
    cor_ersp = cor_ersp(:,:,:,2:end);
    erspdata1 = cor_ersp;
else
    fprintf('No empty subject in cor_ersp\n')
    erspdata1 = cor_ersp;
end

% error corrected trials
if (sum(squeeze(errcor_ersp(:,:,:,1) ~=0),'all') == 0) && (length(errcor_id) ~= size(errcor_ersp,4))
    fprintf('First subject in error corrected ersp data is empty\nDeleting errcor_ersp(:,:,:,1)\n')
    errcor_ersp = errcor_ersp(:,:,:,2:end);
    erspdata2 = errcor_ersp;
else
    fprintf('No empty subject in errcor_ersp\n')
    erspdata2 = errcor_ersp;
end

numChans = size(cor_ersp,1);
numFreqs = size(cor_ersp,2);
numTimes = size(cor_ersp,3);
if exist('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/two_sample_ttest_cor_errcor.mat','file')
    fprintf('already computed two-sample t-test for correct and error corrected trials\n')
else
    % setting up progress bar
    f = waitbar(0,'Chan','Name','Two-sample T-test Channel Progress');
    for i=1:numChans
        % updating waitbar
        waitbar(i/numChans,f,sprintf('Channel %d/%d',i,numChans))
        for j = 1:numFreqs
            for k = 1:numTimes
                data1 = squeeze(erspdata1(i,j,k,:));
                data2 = squeeze(erspdata2(i,j,k,:));
                [~,p(i,j,k),~,stats] = ttest2(data1,data2);
                t(i,j,k) = stats.tstat;
            end
        end    
    end
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/two_sample_ttest_cor_errcor.mat','p','t')
end