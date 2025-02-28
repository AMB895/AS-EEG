%% Plot stats on ERSP data function
function [] = plot_ersp_stats(channels,test,trialtype,mask)
% inputs:
%   channels- vector of channel numbers to plot
%   test- string of the stats test you want to plot:
%           'ttest' = one sample t-test
%           'invage' = linear regression with inverse age
%           'age' = linear regression with age
%   trialtype - 0,1,or 2
%           0 = incorrect vs. correct trials
%           1 = correct trials
%           2 = error corrected vs. correct trials
%   mask - 0 or 1
%           0 = plot with no significance mask
%           1 = plot with significance mask

% Determining number of channels to plot
numChansToPlot = length(channels);

% Getting times and frequencies
% all ERSP correct trials path
allerspcorrectdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/CorrectTrials';
allerspcorrectdataname = 'all_ERSP_DATA_cor.mat';
ALLDATA = load(fullfile(allerspcorrectdatapath,allerspcorrectdataname)); % for times and freqs
times = ALLDATA.times;
freqs = ALLDATA.freqs;

if trialtype == 1
    % all ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/CorrectTrials';
    % Files being loaded
    ttestdataname = 'ttest_outputs_cor.mat';
    ageregressdataname = 'linregress_age_outputs_cor.mat';
    invageregressdataname = 'linregress_invage_outputs_cor.mat';
    subtitstr = 'Correct Trials';
elseif trialtype == 0
    % all ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/';
    % Files being loaded
    ttestdataname = 'two_sample_ttest_cor_incor.mat';
    subtitstr = 'Incorrect vs. Correct Trials';
    if strcmp(test,'age') || strcmp(test,'invage')
        error("Did not compute age or 1/age effects for correct vs. incorrect trials. Use 'ttest'\n")
    end
elseif trialtype == 2
    % all ERSP path
    allerspdatapath = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodPowerAnalysis/';
    % Files being loaded
    ttestdataname = 'two_sample_ttest_cor_errcor.mat';
    subtitstr = 'Error Corrected vs. Correct Trials';
    if strcmp(test,'age') || strcmp(test,'invage')
        error(sprintf('Did not compute age or 1/age effects for correct vs. error corrected trials.\nUse "ttest" as test to plot\n'))
    end
end

% determine type of test to plot
if strcmp(test, 'ttest')
    % loading t-test data
    try
        ALLTTESTDATA = load(fullfile(allerspdatapath,ttestdataname)); % contains p vals and t-statistics
        if trialtype == 1
            tval = ALLTTESTDATA.tval;
            pval_ttest = ALLTTESTDATA.pval_ttest;
        elseif trialtype == 0 || trialtype == 2
            tval = ALLTTESTDATA.t;
            pval_ttest = ALLTTESTDATA.p;
        end
    catch e
        error('Run ersp_stats.m to get t-test outputs')
    end
    
    for currentChan = 1:numChansToPlot
        % defining current channel data
        chanName = string(channels(currentChan));
        figtitle = strcat("Channel ",chanName," Group Activation");
        chan_t = squeeze(tval(channels(currentChan),:,:));
        
        if mask
            % plotting with significance mask
            chan_p = squeeze(pval_ttest(channels(currentChan),:,:));
            sigmask = chan_p < 0.05;
            figure;
            surf(times,freqs,sigmask.*chan_t,'EdgeColor','none')
            view(2); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
            title(figtitle); subtitle(subtitstr)
            colormap('parula'); clim([-10 10]); colorbar
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
        else
            % plotting without significance mask
            figure;
            surf(times,freqs,chan_t,'EdgeColor','none')
            view(2); xlabel('Time (ms)'); ylabel('Frequency (Hz)');
            title(figtitle); subtitle(subtitstr)
            colormap('parula'); clim([-10 10]); colorbar
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
        end
    end
    
elseif strcmp(test, 'age')
    % loading linear regression with age data
    try
        ALLAGEDATA = load(fullfile(allerspdatapath,ageregressdataname)); % contains coefficients, p vals and F-statistics
        b_age = ALLAGEDATA.b_age;
        fval = ALLAGEDATA.fval;
        pval = ALLAGEDATA.pval;
    catch e
        error('Run ersp_stats.m to get linear regression with age outputs')
    end
    
    for currentChan = 1:numChansToPlot
        % defining current channel data
        chanName = string(channels(currentChan));
        figtitle = strcat("Channel ",chanName," Age Effects");
        chan_b = squeeze(b_age(channels(currentChan),:,:));
        chan_F = squeeze(fval(channels(currentChan),:,:));
        
        if mask
            % plotting with mask
            chan_p = squeeze(pval(channels(currentChan),:,:));
            sigmask = chan_p < 0.05;
            
            figure;
            subplot(1,2,1)
            surf(times,freqs,sigmask.*chan_b,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([-0.07 0.07]); colorbar
            title('Age Coefficient'); subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            subplot(1,2,2)
            surf(times,freqs,sigmask.*chan_F,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([0 20]); colorbar
            title('F-statistic'); subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            sgtitle(figtitle)
        else
             % plotting without mask
            figure;
            subplot(1,2,1)
            surf(times,freqs,chan_b,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([-0.07 0.07]); colorbar
            title('Age Coefficient');subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            subplot(1,2,2)
            surf(times,freqs,chan_F,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([0 20]); colorbar
            title('F-statistic');subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            sgtitle(figtitle)
        end    
    end
    
elseif strcmp(test,'invage')
    % loading linear regression with inverse age data
    try
        ALLINVAGEDATA = load(fullfile(allerspdatapath,invageregressdataname)); % contains coefficients, p vals and F-statistics
        b_age_inv = ALLINVAGEDATA.b_age_inv;
        fval_inv = ALLINVAGEDATA.fval_inv;
        pval_inv = ALLINVAGEDATA.pval_inv;
    catch e
        error('Run ersp_stats.m to get linear regression with inverse age outputs')
    end
    
    for currentChan = 1:numChansToPlot
        % defining current channel data
        chanName = string(channels(currentChan));
        figtitle = strcat("Channel ",chanName," 1/Age Effects");
        chan_b = squeeze(b_age_inv(channels(currentChan),:,:));
        chan_b = -1.*chan_b;
        chan_F = squeeze(fval_inv(channels(currentChan),:,:));
        
        if mask
            % plotting with mask
            chan_p = squeeze(pval_inv(channels(currentChan),:,:));
            sigmask = chan_p < 0.05;
            
            figure;
            subplot(1,2,1)
            surf(times,freqs,sigmask.*chan_b,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([-20 20]); colorbar
            title('1/Age Coefficient');subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            subplot(1,2,2)
            surf(times,freqs,sigmask.*chan_F,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([0 20]); colorbar
            title('F-statistic');subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            sgtitle(figtitle)
        else
             % plotting without mask
            figure;
            subplot(1,2,1)
            surf(times,freqs,chan_b,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([-0.07 0.07]); colorbar
            title('Age Coefficient');subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            subplot(1,2,2)
            surf(times,freqs,chan_F,'EdgeColor','none')
            view(2)
            xlabel('Time (ms)'); ylabel('Frequency (Hz)')
            colormap('parula'); clim([0 20]); colorbar
            title('F-statistic');subtitle(subtitstr)
            
            % plotting fixation cross onset and stimulus onset
            hold on
            xline(0,'--r','LineWidth',1.5)
            xline(500,'--y','LineWidth',1.5)
            
            sgtitle(figtitle)
        end    
    end
    
else
    error(sprintf("Incorrect string for test to plot\nString must be one of the following:\n'ttest'\n'age'\n'invage'\n"))
end

end
% %% Plotting results
% % Loading outputs and ERSP data
% load(fullfile(allerspdatapath,'one_sample_ttest_outputs.mat'))
% load(fullfile(allerspdatapath,'inv_age_linear_regress_outputs.mat'))
% load(fullfile(allerspdatapath,'linear_regress_outputs.mat'))
% load(fullfile(allerspdatapath,'ERSP_DATA.mat')) % to get times and frequencies
% %% Group Activation (One-sample T-test)
% numChans = size(tval,1);
%  
% % where to save figures
% groupactpath = fullfile(allerspdatapath,'GroupActivationPlots/');
% 
% % plotting t-test results for each channel
% for currentChan = 1:numChans
%     
%     close all
%     % defining current channel data
%     chanName = string(currentChan);
%     figtitle = strcat("Channel ",chanName," Group Activation");
%     figsavename = strcat('Chan',chanName,'.pdf');
%     chan_tval = squeeze(tval(currentChan,:,:));
%     chan_h = squeeze(h(currentChan,:,:));
%     
%     % plotting t val with h mask
%     fig1 = figure;
%     surf(times,freqs,chan_h.*chan_tval,'EdgeColor','none')
%     view(2)
%     title(figtitle); subtitle('t-statistic with significance mask')
%     xlabel('Time (ms)'); ylabel('Frequency (Hz)')
%     colormap('parula'); clim([-10 10]); colorbar
%     
%     % line at onset of red fixation cross
%     hold on
%     xline(0,'--r','LineWidth',1.5)
% 
%     % line at onset of yellow target
%     xline(500,'--y','LineWidth',1.5)
%     
%     % saving figure
%     saveas(fig1,fullfile(groupactpath,figsavename))
%     
% end
% 
% %% Age-effects Linear Regression
% numChans = size(fval,1);
% % where to save figures
% ageeffectspath = fullfile(allerspdatapath,'AgeEffectsPlots/');
% 
% % plotting linear regression results for each channel
% for currentChan = 1:numChans
%     
%     close all
%     % defining current channel data
%     chanName = string(currentChan);
%     figtitle = strcat("Channel ",chanName," Age Effects- Linear Regression");
%     figsavename = strcat('Chan',chanName,'.pdf');
%     % defining current channel data
%     chanName = string(currentChan);
%     figtitle = strcat("Channel ",chanName," Age Effects- Linear Regression");
%     figsavename = strcat('Chan',chanName,'.pdf');
%     chan_bage = squeeze(b_age(currentChan,:,:));
%     chan_pval = squeeze(pval(currentChan,:,:));
%     chan_fval = squeeze(fval(currentChan,:,:));
%     sig_mask = chan_pval < 0.05;
%     
%     % plotting age coefficient (b) with significance mask
%     fig2 = figure;
%     subplot(1,2,1)
%     surf(times,freqs,sig_mask.*chan_bage,'EdgeColor','none')
%     view(2)
%     subtitle('Age Coefficient with significance mask')
%     xlabel('Time (ms)'); ylabel('Frequency (Hz)')
%     colormap('parula'); clim([-0.07 0.07]); colorbar
% 
%     % line at onset of red fixation cross
%     hold on
%     xline(0,'--r','LineWidth',1.5)
%     
%     % line at onset of yellow target
%     xline(500,'--y','LineWidth',1.5)
% 
%     % plotting f values with significance mask
%     subplot(1,2,2)
%     surf(times,freqs,sig_mask.*chan_fval,'EdgeColor','none')
%     view(2)
%     subtitle('F-statistic with significance mask')
%     xlabel('Time (ms)'); ylabel('Frequency (Hz)')
%     colormap('parula'); clim([0 20]); colorbar
% 
%      % line at onset of red fixation cross
%     hold on
%     xline(0,'--r','LineWidth',1.5)
%     
%     % line at onset of yellow target
%     xline(500,'--y','LineWidth',1.5)
% 
%     sgtitle(figtitle)
%     
%     % saving figure
%     fig2.PaperOrientation = 'landscape';
%     fig2.PaperPosition = [0.25 3 10 5];
%     saveas(fig2,fullfile(ageeffectspath,figsavename))
%     
% end
% 
% %% Inverse Age Effects - Linear Regression
% numChans = size(fval_inv,1);
% % where to save figures
% invageeffectspath = fullfile(allerspdatapath,'InverseAgeEffectsPlots/');
% 
% % plotting linear regression results for each channel
% for currentChan = 1:numChans
%     
%     close all
%     % defining current channel data
%     chanName = string(currentChan);
%     figtitle = strcat("Channel ",chanName," 1/Age Effects- Linear Regression");
%     figsavename = strcat('Chan',chanName,'.pdf');
%     chan_bageinv = squeeze(b_age_inv(currentChan,:,:));
%     fix_bageinv = -1.*chan_bageinv;
%     chan_pvalinv = squeeze(pval_inv(currentChan,:,:));
%     chan_fvalinv = squeeze(fval_inv(currentChan,:,:));
%     sig_mask = chan_pvalinv < 0.05;
%     
%     % plotting age coefficient (b) with significance mask
%     fig3 = figure;
%     subplot(1,2,1)
%     surf(times,freqs,sig_mask.*fix_bageinv,'EdgeColor','none')
%     view(2)
%     subtitle('1/Age Coefficient with significance mask')
%     xlabel('Time (ms)'); ylabel('Frequency (Hz)')
%     colormap('parula'); clim([-20 20]); colorbar
%     
%     % line at onset of red fixation cross
%     hold on
%     xline(0,'--r','LineWidth',1.5)
%     
%     % line at onset of yellow target
%     xline(500,'--y','LineWidth',1.5)
%     
%     % plotting f values with significance mask
%     subplot(1,2,2)
%     surf(times,freqs,sig_mask.*chan_fvalinv,'EdgeColor','none')
%     view(2)
%     subtitle('F-statistic with significance mask')
%     xlabel('Time (ms)'); ylabel('Frequency (Hz)')
%     colormap('parula'); clim([0 20]); colorbar
% 
%     % line at onset of red fixation cross
%     hold on
%     xline(0,'--r','LineWidth',1.5)
%     
%     % line at onset of yellow target
%     xline(500,'--y','LineWidth',1.5)
% 
%     sgtitle(figtitle)
%     
%     % saving figure
%     fig3.PaperOrientation = 'landscape';
%     fig3.PaperPosition = [0.25 3 10 5];
%     saveas(fig3,fullfile(invageeffectspath,figsavename))
%     
% end
