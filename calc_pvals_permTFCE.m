%% Calculating p-values from Permutation TFCE analysis
% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/TFCE/')
% load original TFCE scores
load('TFCE_scores_correct.mat')
load('TFCE_scores_cor_errcor.mat')
load('TFCE_scores_cor_errcor.mat')

% load maximum permuted TFCE scores
load('permutationTFCE_correct.mat')
load('permutationTFCE_cor_incor.mat')
load('permutationTFCE_cor_errcor.mat')

% plot null distributions for each channel
% plotdata = 1; % plot
plotdata = 0; % do not plot
%% Method 1: (per channel)
%   1. Count number of iterations that produce permuted TFCE > max original TCFE
%   2. Divide count by total number of permutations = p-value
nChans = size(tfce_scores_group_cor,1);

for currentChan = 1:nChans
    % Correct Trials Group Activation
    maxOGtfce_cor_group = max(tfce_scores_group_cor(currentChan,:,:),[],'all');
    permtfce_cor_group = cor_max_perm_tfce_group(currentChan,:,:);
    p_cor_group(currentChan) = nnz(permtfce_cor_group > maxOGtfce_cor_group)/1000;
    % Correct Trials Age Effects
    maxOGtfce_cor_age = max(tfce_scores_age_cor(currentChan,:,:),[],'all');
    permtfce_cor_age = cor_max_perm_tfce_age(currentChan,:,:);
    p_cor_age(currentChan) = nnz(permtfce_cor_age > maxOGtfce_cor_age)/1000;
    % Correct Trials Inverse Age Effects
    maxOGtfce_cor_ageinv = max(tfce_scores_ageinv_cor(currentChan,:,:),[],'all');
    permtfce_cor_ageinv = cor_max_perm_tfce_ageinv(currentChan,:,:);
    p_cor_ageinv(currentChan) = nnz(permtfce_cor_ageinv > maxOGtfce_cor_ageinv)/1000;
    % Correct vs. Incorrect Group Activation
    maxOGfce_cor_incor_group = max(tfce_scores_cor_incor(currentChan,:,:),[],'all');
    permtfce_cor_incor_group = cor_incor_max_perm_tfce_group(currentChan,:,:);
    p_cor_incor_group(currentChan) = nnz(permtfce_cor_incor_group > maxOGfce_cor_incor_group)/1000;
     % Correct vs. Error Corrected Group Activation
    maxOGfce_cor_errcor_group = max(tfce_scores_cor_errcor(currentChan,:,:),[],'all');
    permtfce_cor_errcor_group = cor_errcor_max_perm_tfce_group(currentChan,:,:);
    p_cor_errcor_group(currentChan) = nnz(permtfce_cor_errcor_group > maxOGfce_cor_errcor_group)/1000;
end

%% Method 2: (per Channel)
% 1. Generate null distribution of max permuted TFCE scores
% 2. Find 95th percentile
% 3. Threshold original TFCE scores at 95th threshold

nChans = size(tfce_scores_group_cor,1);
for currentChan = 1:nChans
    close all
    chanString = string(currentChan);
    % defining current channel's maximum permuted tfce score
    chan_max_perm_cor_group = cor_max_perm_tfce_group(currentChan,:,:);
    chan_max_perm_cor_age = cor_max_perm_tfce_age(currentChan,:,:);
    chan_max_perm_cor_ageinv = cor_max_perm_tfce_ageinv(currentChan,:,:);
    chan_max_perm_cor_incor_group = cor_incor_max_perm_tfce_group(currentChan,:,:);
    chan_max_perm_cor_errcor_group = cor_errcor_max_perm_tfce_group(currentChan,:,:);

    % determining 95th percentile
    n_corGroup = length(chan_max_perm_cor_group)*0.95;
    n_corAge = length(chan_max_perm_cor_age)*0.95;
    n_corAgeInv = length(chan_max_perm_cor_ageinv)*0.95;
    n_cor_incorGroup = length(chan_max_perm_cor_incor_group)*0.95;
    n_cor_errcorGroup = length(chan_max_perm_cor_errcor_group)*0.95;

    % sort maximum permuted tfce scores in ascending order
    sorted_corGroup = sort(chan_max_perm_cor_group);
    sorted_corAge = sort(chan_max_perm_cor_age);
    sorted_corAgeInv = sort(chan_max_perm_cor_ageinv);
    sorted_cor_incorGroup = sort(chan_max_perm_cor_incor_group);
    sorted_cor_errcorGroup = sort(chan_max_perm_cor_errcor_group);

    % determine 95th percentile threshold
    thres95_corGroup(currentChan) = chan_max_perm_cor_group(n_corGroup);
    thres95_corAge(currentChan) = chan_max_perm_cor_age(n_corAge);
    thres95_corAgeInv(currentChan) = chan_max_perm_cor_ageinv(n_corAgeInv);
    thres95_cor_incorGroup(currentChan) = chan_max_perm_cor_incor_group(n_cor_incorGroup);
    thres95_cor_errcorGroup(currentChan) = chan_max_perm_cor_errcor_group(n_cor_errcorGroup);

    % threshold mask
    mask_corGroup(currentChan) = chan_max_perm_cor_group > thres95_corGroup(currentChan);
    mask_corAge(currentChan) = chan_max_perm_cor_age > thres95_corAge(currentChan);
    mask_corAgeInv(currentChan) = chan_max_perm_cor_ageinv > thres95_corAgeInv(currentChan);
    mask_cor_incorGroup(currentChan) = chan_max_perm_cor_incor_group > thres95_cor_incorGroup(currentChan);
    mask_cor_errcorGroup(currentChan) = chan_max_perm_cor_errcor_group > thres95_cor_errcorGroup(currentChan);
    
    if plotdata
        % plot thresholded TFCE plots for each channel
        % Correct Trials - Group Activation
        figure;
        surf(times,freqs,mask_corGroup.*squeeze(tfce_scores_group_cor(currentChan,:,:)),'EdgeColor','none');
        view(2); colorbar; 
        hold on
        xline(0,'--r','LineWidth',1.5); xline(500,'--y','LineWidth',1.5);
        xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
        title('Thresholded TFCE Correct Trials - Group Activation'); subtitle(sprintf('Channel %s',chanString));
    
        % Correct Trials - Age Effects
        figure;
        surf(times,freqs,mask_corAge.*squeeze(tfce_scores_age_cor(currentChan,:,:)),'EdgeColor','none');
        view(2); colorbar;
        hold on
        xline(0,'--r','LineWidth',1.5); xline(500,'--y','LineWidth',1.5);
        xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
        title('Thresholded TFCE Correct Trials - Age Effects'); subtitle(sprintf('Channel %s',chanString));
    
    
        % Correct Trials - Inverse Age Effects
        figure;
        surf(times,freqs,mask_corAgeInv.*squeeze(tfce_scores_ageinv_cor(currentChan,:,:)),'EdgeColor','none');
        view(2); colorbar;
        hold on
        xline(0,'--r','LineWidth',1.5); xline(500,'--y','LineWidth',1.5);
        xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
        title('Thresholded TFCE Correct Trials - Inverse Age Effects'); subtitle(sprintf('Channel %s',chanString));
        
        % Correct vs. Incorrect Trials - Group Activation
        figure;
        surf(times,freqs,mask_cor_incorGroup.*squeeze(tfce_scores_group_cor_incor(currentChan,:,:)),'EdgeColor','none');
        view(2); colorbar;
        hold on
        xline(0,'--r','LineWidth',1.5); xline(500,'--y','LineWidth',1.5);
        xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
        title('Thresholded TFCE Correct vs. Incorrect Trials - Group Activation'); subtitle(sprintf('Channel %s',chanString));
        
    
        % Correct vs. Error Corrected Trials - Group Activation
        figure;
        surf(times,freqs,mask_cor_errcorGroup.*squeeze(tfce_scores_group_cor_errcor(currentChan,:,:)),'EdgeColor','none');
        view(2); colorbar;
        hold on
        xline(0,'--r','LineWidth',1.5); xline(500,'--y','LineWidth',1.5);
        xlabel('Time (ms)'); ylabel('Frequency (Hz)'); 
        title('Thresholded TFCE Correct vs. Error Corrected Trials - Group Activation'); subtitle(sprintf('Channel %s',chanString));
    
        % plotting null distribution of maximum permuted TFCE scores
        % Correct Trials - Group Activation
        figure;
        histogram(chan_max_perm_cor_group)
        hold on
        xline(thres95_corGroup(currentChan),'-r','LineWidth',1.5)
        xlabel('Max Permuted TFCE Score'); ylabel('Count');
        title('Null Distribution Correct Trials - Group Activation'); subtitle(sprintf('Channel %s',chanString));
        
        % Correct Trials - Age Effects
        figure;
        histogram(chan_max_perm_cor_age)
        hold on
        xline(thres95_corAge(currentChan),'-r','LineWidth',1.5)
        xlabel('Max Permuted TFCE Score'); ylabel('Count');
        title('Null Distribution Correct Trials - Age Effects'); subtitle(sprintf('Channel %s',chanString));
    
        % Correct Trials - Inverse Age Effects
        figure;
        histogram(chan_max_perm_cor_ageinv)
        hold on
        xline(thres95_corAgeInv(currentChan),'-r','LineWidth',1.5)
        xlabel('Max Permuted TFCE Score'); ylabel('Count');
        title('Null Distribution Correct Trials - Inverse Age Effects'); subtitle(sprintf('Channel %s',chanString));
    
        % Correct vs. Incorrect Trials - Group Activation
        figure;
        histogram(chan_max_perm_cor_incor_group)
        hold on
        xline(thres95_cor_incorGroup(currentChan),'-r','LineWidth',1.5)
        xlabel('Max Permuted TFCE Score'); ylabel('Count');
        title('Null Distribution Correct vs. Incorrect Trials - Group Activation'); subtitle(sprintf('Channel %s',chanString));
    
        % Correct vs. Error Corrected Trials - Group Activation
        figure;
        histogram(chan_max_perm_cor_errcor_group)
        hold on
        xline(thres95_cor_errcorGroup(currentChan),'-r','LineWidth',1.5)
        xlabel('Max Permuted TFCE Score'); ylabel('Count');
        title('Null Distribution Correct vs. Error Corrected Trials - Group Activation'); subtitle(sprintf('Channel %s',chanString));
    else
        fprintf('Not plotting thresholded TFCE plots; set plotdata=1 to plot\n')
    end
end

% saving p-values calculated in method 1, TFCE threshold score at 0.05
% significance level and the threshold masks at 0.05 significance
save('/Volumes/Hera/Abby/AS_EEG/TFCE/TFCEpvalues.mat','p_cor_errcor_group','p_cor_incor_group','p_cor_ageinv','p_cor_age','p_cor_group')
save('/Volumes/Hera/Abby/AS_EEG/TFCE/thresholds95.mat','thres95_cor_errcorGroup','thres95_cor_incorGroup','thres95_corAgeInv','thres95_corAge','thres95_corGroup')
save('/Volumes/Hera/Abby/AS_EEG/TFCE/thres95masks.mat','mask_corGroup','mask_corAgeInv','mask_corAge','mask_cor_errcorGroup','mask_cor_incorGroup')