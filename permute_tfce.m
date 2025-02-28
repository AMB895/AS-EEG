%% Permutation Test to get p-values for TFCE scores
% 1. Permute t-values/ F-values 
% 2. Run TFCE on permuted values
% 3. Find max TFCE score for each iteration
%% DO THIS ON SUPER COMPUTER!!!
%% Setting up
% increase number of cores for parallel processing
% limo defaults to 35, but only see 22
setenv('NUMBER_OF_PROCESSORS','50')
% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')

% start eeglab
[ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;

% TFCE data path
tfcepath = '/Volumes/Hera/Abby/AS_EEG/TFCE/';

% Names and paths of data to be loaded
maindir = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/';
cortrialspath = [maindir,'CorrectTrials/'];
corttestname = 'ttest_outputs_cor.mat';
coragename = 'linregress_age_outputs_cor.mat';
corageinvname = 'lineregress_invage_outputs_cor.mat';
corincorttestname = 'two_sample_ttest_cor_incor.mat';
corerrcorttestname = 'two_sample_ttest_cor_errcor.mat';

% names of files to be created
cor_permtfcename = 'permutationTFCE_correct.mat';
cor_incor_permtfcename = 'permutationTFCE_cor_incor.mat';
cor_errcor_permtfcename = 'permutationTFCE_cor_errcor.mat';

nperm = 1000;
f = waitbar(0,'Starting','Name','TFCE Permutation Progress');

% load original t/F values: channels x freqs x times
COR_TVAL = load(fullfile(cortrialspath,corttestname),'tval');
cor_tval = COR_TVAL.tval;
COR_FVAL_AGE = load(fullfile(cortrialspath,coragename),'fval');
cor_fval_age = COR_FVAL_AGE.fval;
COR_FVAL_AGEINV = load(fullfile(cortrialspath,corageinvname),'fval_inv');
cor_fval_ageinv = COR_FVAL_AGEINV.fval_inv;

COR_INCOR_TTEST = load(fullfile(maindir,corincorttestname),'t');
cor_incor_tval = COR_INCOR_TTEST.t;
COR_ERRCOR_TTEST = load(fullfile(maindir,corerrcorttestname),'t');
cor_errcor_tval = COR_ERRCOR_TTEST.t;

% Get number of channels to loop through
nchans = size(cor_tval,1);

% preallocating max permuted tfce score vectors for speed
cor_max_perm_tfce_group = zeros(nchans,nperm);
cor_max_perm_tfce_age = zeros(nchans,nperm);
cor_max_perm_tfce_ageinv = zeros(nchans,nperm);

cor_incor_max_perm_tfce_group = zeros(nchans,nperm);
cor_errcor_max_perm_tfce_group = zeros(nchans,nperm);

for n=1:nperm
    %% Correct Trials
    % getting total size of t-values matrix
    cor_totalSize = size(cor_tval,1)*size(cor_tval,2)*size(cor_tval,3);
    
    % reshaping tfce scores to a 1 x totalSize vector
    cor_reshaped_t = reshape(cor_tval,[1,cor_totalSize]);
    cor_reshaped_f = reshape(cor_fval_age,[1, cor_totalSize]);
    cor_reshaped_finv = reshape(cor_fval_ageinv,[1,cor_totalSize]);
    
    % Permute t-value matrix
    cor_perm_t = cor_reshaped_t(randperm(length(cor_reshaped_t)));
    cor_perm_f = cor_reshaped_f(randperm(length(cor_reshaped_f)));
    cor_perm_finv = cor_reshaped_finv(randperm(length(cor_reshaped_finv)));
    
    % Reshape permuted t-values back to original size
    cor_perm_t = reshape(cor_perm_t,[size(cor_tval,1) size(cor_tval,2) size(cor_tval,3)]);
    cor_perm_f = reshape(cor_perm_f,[size(cor_tval,1) size(cor_tval,2) size(cor_tval,3)]);
    cor_perm_finv = reshape(cor_perm_finv, [size(cor_tval,1) size(cor_tval,2) size(cor_tval,3)]);
    
    % updating waitbar
    waitbar(n/nperm,f,sprintf('%d',n))
    
    % Run TFCE on permuted t-values/ F-values for correct trials
    % limo_tfce inputs:
    %   type: 3 (for 3D data) where data is channels x freqs x times
    %   data: 3D data
    %   channeighbmatrix: logical matrix indicating neighboring channels
    %   updatebar: 0 or 1 to display TFCE thresholding progress
    %               0 = do not show progress bar
    %               1 = show progress bar
    cor_perm_tfce_score_group = limo_tfce(3,cor_perm_t,channeighbmatrix,0);
    cor_perm_tfce_score_age = limo_tfce(3,cor_perm_f,channeighbmatrix,0);
    cor_perm_tfce_score_ageinv = limo_tfce(3,cor_perm_finv,channeighbmatrix,0);

    % Save maximum TFCE score over all channels
    cor_max_perm_tfce_group(:,n) = max(cor_perm_tfce_score_group,[],[2 3]);
    cor_max_perm_tfce_age(:,n) = max(cor_perm_tfce_score_age,[],[2,3]);
    cor_max_perm_tfce_ageinv(:,n) = max(cor_perm_tfce_score_ageinv,[],[2 3]);

    %% Correct and Incorrect Trials
    % getting total size of t-values matrix
    cor_incor_totalSize = size(cor_incor_tval,1)*size(cor_incor_tval,2)*size(cor_incor_tval,3);
    
    % reshaping tfce scores to a 1 x totalSize vector
    cor_incor_reshaped_t = reshape(cor_incor_tval,[1,cor_incor_totalSize]);
    
    % Permute t-value matrix
    cor_incor_perm_t = cor_incor_reshaped_t(randperm(length(cor_incor_reshaped_t)));
    
    % Reshape permuted t-values back to original size
    cor_incor_perm_t = reshape(cor_incor_perm_t,[size(cor_incor_tval,1) size(cor_incor_tval,2) size(cor_incor_tval,3)]);
    
    % Run TFCE on permuted t-values/ F-values for correct trials
    cor_incor_perm_tfce_score_group = limo_tfce(3,cor_incor_perm_t,channeighbmatrix,0);

    % Save maximum TFCE score over all channels
    cor_incor_max_perm_tfce_group(:,n) = max(cor_incor_perm_tfce_score_group,[],[2 3]);

    %% Correct and Error Corrected Trials
    % getting total size of t-values matrix
    cor_errcor_totalSize = size(cor_errcor_tval,1)*size(cor_errcor_tval,2)*size(cor_errcor_tval,3);
    
    % reshaping tfce scores to a 1 x totalSize vector
    cor_errcor_reshaped_t = reshape(cor_errcor_tval,[1,cor_errcor_totalSize]);
    
    % Permute t-value matrix
    cor_errcor_perm_t = cor_errcor_reshaped_t(randperm(length(cor_errcor_reshaped_t)));
    
    % Reshape permuted t-values back to original size
    cor_errcor_perm_t = reshape(cor_errcor_perm_t,[size(cor_errcor_tval,1) size(cor_errcor_tval,2) size(cor_errcor_tval,3)]);
    
    % Run TFCE on permuted t-values/ F-values for correct trials
    cor_errcor_perm_tfce_score_group = limo_tfce(3,cor_errcor_perm_t,channeighbmatrix,0);

    % Save maximum TFCE score over all channels
    cor_errcor_max_perm_tfce_group(:,n) = max(cor_errcor_perm_tfce_score_group,[],[2 3]);

end
save(fullfile(tfcepath,cor_permtfcename),'cor_max_perm_tfce_group','cor_max_perm_tfce_age','cor_max_perm_tfce_ageinv')
save(fullfile(tfcepath,cor_incor_permtfcename),'cor_incor_max_perm_tfce_group')
save(fullfile(tfcepath,cor_errcor_permtfcename),'cor_errcor_max_perm_tfce_group')

% save(fullfile(tfcepath,'permutationTFCE.mat'),'max_perm_tfce_groupact','max_perm_tfce_age','max_perm_tfce_ageinv')
% 
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







