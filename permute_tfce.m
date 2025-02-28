%% Permutation Test to get p-values for TFCE scores
% 1. Permute t-values/ F-values 
% 2. Run TFCE on permuted values
% 3. Find max TFCE score for each iteration
%% Setting up
% get host name
[~,hostname] = system('hostname');

% directories and paths for running on rhea
if strcmp(hostname,'rhea.wpic.upmc.edu')
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
    
    % setting up waitbar- only when running on rhea
    f = waitbar(0,'Starting','Name','TFCE Permutation Progress');

elseif strcmp(hostname,'br014.ib.bridges2.psc.edu') % directories & paths for PSC
    % increase number of cores for parallel processing
    % limo defaults to 35, but only see 22
    setenv('NUMBER_OF_PROCESSORS','50')
    % adding necessary paths 
    % FIXME!! rsync toolboxes to PSC
    addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/tools/')
    
    % start eeglab
    [ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;
    
    % TFCE data path
    tfcepath = '/ocean/projects/soc230004p/shared/antisaccade_eeg/';
    
    % Names and paths of data to be loaded
    maindir = '/ocean/projects/soc230004p/shared/antisaccade_eeg/data/';
    corttestname = 'ttest_outputs_cor.mat';
    coragename = 'linregress_age_outputs_cor.mat';
    corageinvname = 'lineregress_invage_outputs_cor.mat';
    corincorttestname = 'two_sample_ttest_cor_incor.mat';
    corerrcorttestname = 'two_sample_ttest_cor_errcor.mat';
    
    % names of files to be created
    cor_permtfcename = 'permutationTFCE_correct.mat';
    cor_incor_permtfcename = 'permutationTFCE_cor_incor.mat';
    cor_errcor_permtfcename = 'permutationTFCE_cor_errcor.mat';
    
    % load original t/F values: channels x freqs x times
    COR_TVAL = load(fullfile(maindir,corttestname),'tval');
    cor_tval = COR_TVAL.tval;
    COR_FVAL_AGE = load(fullfile(maindir,coragename),'fval');
    cor_fval_age = COR_FVAL_AGE.fval;
    COR_FVAL_AGEINV = load(fullfile(maindir,corageinvname),'fval_inv');
    cor_fval_ageinv = COR_FVAL_AGEINV.fval_inv;
    
    COR_INCOR_TTEST = load(fullfile(maindir,corincorttestname),'t');
    cor_incor_tval = COR_INCOR_TTEST.t;
    COR_ERRCOR_TTEST = load(fullfile(maindir,corerrcorttestname),'t');
    cor_errcor_tval = COR_ERRCOR_TTEST.t;

elseif strcmp(hostname,'oaclf1lunalin1.acct.upmchs.net') % cannot run on local computer
    error('Run permute_tfce.m on supercomputer or rhea')
end

% Number of permutations
nperm = 1000;

% Get number of channels
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
    if strcmp(hostname,'rhea.wpic.upmc.edu')
        waitbar(n/nperm,f,sprintf('%d',n))
    end
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

