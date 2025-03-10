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
    mask_corGroup = chan_max_perm_cor_group > thres95_corGroup(currentChan);
    mask_corAge = chan_max_perm_cor_age > thres95_corAge(currentChan);
    mask_corAgeInv = chan_max_perm_cor_ageinv > thres95_corAgeInv(currentChan);
    mask_cor_incorGroup = chan_max_perm_cor_incor_group > thres95_cor_incorGroup(currentChan);
    mask_cor_errcorGroup = chan_max_perm_cor_errcor_group > thres95_cor_errcorGroup(currentChan);

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
end

%%
% %% Visualizing Null Distribution of Permuted TFCE Scores
% % Method 1 to get p-values: 
% %   1. Count number of iterations that produce TFCE > max TFCE
% %   2. Divide count by total number of permutations = p-value
% 
% % Method 2: 
% %   1. generate null distribution of permuted max tfce scores
% %   2. find 95th percentile
% %   3. determine if original tfce scores are greater than 95th percentile threshold
% 
% % load permutation TFCE scores
% load(fullfile(tfcepath,permtfcename))
% 
% nChans = size(max_perm_tfce_ageinv,1);
% 
% for currentChan = 1:nChans
%     % defining channel number for chart title
%     subtitstr = sprintf('Channel %d',currentChan);
% 
%     % defining current channels permutation scores
%     chan_permTFCE_group = max_perm_tfce_groupact(currentChan,:);
%     chan_permTFCE_age = max_perm_tfce_age(currentChan,:);
%     chan_permTFCE_ageinv = max_perm_tfce_ageinv(currentChan,:);
% 
%     % getting current channel's original TFCE scores
%     chan_tfce_group = squeeze(tfce_scores_groupact(currentChan,:,:));
%     chan_tfce_age = squeeze(tfce_scores_age(currentChan,:,:));
%     chan_tfce_ageinv = squeeze(tfce_scores_ageinv(currentChan,:,:));
% 
%     % finding maximum TFCE score in original data
%     chan_max_tfce_group = max(chan_tfce_group,[],'all');
%     chan_max_tfce_age = max(chan_tfce_age,[],'all');
%     chan_max_tfce_ageinv = max(chan_tfce_ageinv,[],'all');
% 
%     % Method 1: count number of times permuted TFCE score > max TFCE score
%     p_group(currentChan) = length(find(chan_permTFCE_group >= chan_max_tfce_group))/nperm;
%     p_age(currentChan) = length(find(chan_permTFCE_age >= chan_max_tfce_age))/nperm;
%     p_ageinv(currentChan) = length(find(chan_permTFCE_ageinv >= chan_max_tfce_ageinv))/nperm;
% 
%     % determine 95th percentile of permuted tfce scores
%     n_group = length(chan_permTFCE_group)*0.95;
%     n_age =  length(chan_permTFCE_age)*0.95;
%     n_ageinv = length(chan_permTFCE_ageinv)*0.95;
% 
%     % sort max permuted tfce scores in ascending order
%     sorted_tfce_group = sort(chan_permTFCE_group);
%     sorted_tfce_age = sort(chan_permTFCE_age);
%     sorted_tfce_ageinv = sort(chan_permTFCE_ageinv);
% 
%     % determine 95th percentile threshold TFCE value
%     thres95_group(currentChan) = sorted_tfce_group(n_group);
%     thres95_age(currentChan) = sorted_tfce_age(n_age);
%     thres95_ageinv(currentChan) = sorted_tfce_ageinv(n_ageinv);
% 
%     % Method 2: looking at thresholded original TFCE plots
%     mask_group = chan_tfce_group > thres95_group(currentChan);
%     mask_age = chan_tfce_age > thres95_age(currentChan);
%     mask_ageinv = chan_tfce_ageinv > thres95_ageinv(currentChan);
% 
%     % % group activation
%     % figure;
%     % surf(times,freqs,mask_group.*chan_tfce_group,'EdgeColor','none')
%     % view(2); xlabel('Time (ms)'); ylabel('Frequency (Hz)'); colorbar; clim([0 5000])
%     % title('p < 0.05 Group Activation'); subtitle(subtitstr); colormap('turbo')
%     % hold on
%     % xline(0,'--r','LineWidth',1.5)
%     % xline(500,'--y','LineWidth',1.5)
%     % 
%     %  % age effects
%     % figure;
%     % surf(times,freqs,mask_age.*chan_tfce_age,'EdgeColor','none')
%     % view(2); xlabel('Time (ms)'); ylabel('Frequency (Hz)'); colorbar;clim([0 5000])
%     % title('p < 0.05 Age Effects'); subtitle(subtitstr); colormap("turbo")
%     % hold on
%     % xline(0,'--r','LineWidth',1.5)
%     % xline(500,'--y','LineWidth',1.5)
%     % 
%     % % inverse age effects
%     % figure;
%     % surf(times,freqs,mask_ageinv.*chan_tfce_ageinv,'EdgeColor','none')
%     % view(2); xlabel('Time (ms)'); ylabel('Frequency (Hz)'); colorbar;clim([0 5000])
%     % title('p < 0.05 Inverse Age Effects'); subtitle(subtitstr); colormap('turbo')
%     % hold on
%     % xline(0,'--r','LineWidth',1.5)
%     % xline(500,'--y','LineWidth',1.5)
%     % 
%     % % plotting histograms of null distribution from permutations
%     % % group activation
%     % figure;
%     % histogram(chan_permTFCE_group)
%     % xlabel('Max TFCE Score'); ylabel('Count'); title('Distribution of Permuted TFCE Scores: Group Activation'); subtitle(subtitstr);
%     % hold on
%     % xline(thres95_group(currentChan),'r')
%     % 
%     % % % age effects 
%     % figure;
%     % histogram(chan_permTFCE_age)
%     % xlabel('Max TFCE Score'); ylabel('Count'); title('Distribution of Permuted TFCE Scores: Age Effects');subtitle(subtitstr);
%     % hold on
%     % xline(thres95_age(currentChan),'r')
%     % 
%     % % % inverse age effects 
%     % figure;
%     % histogram(chan_permTFCE_ageinv)
%     % xlabel('Max TFCE Score'); ylabel('Count'); title('Distribution of Permuted TFCE Scores: Inverse Age Effects');subtitle(subtitstr);
%     % hold on
%     % xline(thres95_ageinv(currentChan),'r')
%     % 
%     % close all
% end