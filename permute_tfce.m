%% Permutation Test to get p-values for TFCE scores
% 1. Permute t-values/ F-values 
% 2. Run TFCE on permuted values
% 3. Find max TFCE score for each iteration
%% Setting up
% get host name
[~,hostname] = system('hostname');

% directories and paths for running on rhea
if strcmp(hostname(1:end-1),'rhea.wpic.upmc.edu')
    % adding necessary paths
    addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
    addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
    addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
    limo_check_ppool_psc;
    
    % increase number of cores for parallel processing
    setenv('NUMBER_OF_PROCESSORS','1')
    N = 20;
    
    % start eeglab
    [ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;
    
    % TFCE data save path
    tfcepath = '/Volumes/Hera/Abby/AS_EEG/TFCE/';
    
    % Names and paths of data to be loaded
    maindir = '/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/';
    cortrialsdir = [maindir,'CorrectTrials/'];

    % setting up waitbar- only when running on rhea
%     f = waitbar(0,'Starting','Name','TFCE Permutation Progress');

    % load channel neighborhood matrix
    chanStructure = load('/Volumes/Hera/Abby/AS_EEG/STUDY/logicalchanneighbmatrix.mat');
    channeighbmatrix = chanStructure.channeighbmatrix;
    
elseif exist('/ocean/projects/','dir')
   
    % adding necessary paths 
    addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/tools/')

    limo_check_ppool_psc;

    % send NUMBER_OF_PROCESSORS to .bash script
    getenv('NUMBER_OF_PROCESSORS')

    N = 128;
    
    % start eeglab
    [ALLEEG,EEG,CURRENTSET,ALLCOM]=eeglab;
    
    % TFCE data path
    tfcepath = '/ocean/projects/soc230004p/shared/antisaccade_eeg/';
    
    % Names and paths of data to be loaded
    maindir = '/ocean/projects/soc230004p/shared/antisaccade_eeg/data/';
    cortrialsdir = '/ocean/projects/soc230004p/shared/antisaccade_eeg/data/';
 
    % load channel neighborhood matrix
    chanStructure = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/logicalchanneighbmatrix.mat');
    channeighbmatrix = chanStructure.channeighbmatrix;
    
elseif strcmp(hostname(1:end-1),'oaclf1lunalin1.acct.upmchs.net') % cannot run on local computer
    error('Run permute_tfce.m on supercomputer or rhea')
end

% Load t/F values [chans, freqs, times]
COR_TTEST = load(fullfile(cortrialsdir,'ttest_outputs_cor.mat'),'tval');
cor_tval = COR_TTEST.tval;
% COR_AGE = load(fullfile(cortrialsdir,'linregress_age_outputs_cor.mat'),'fval');
% cor_fval_age = COR_AGE.fval;
COR_AGEINV = load(fullfile(cortrialsdir,'linregress_invage_outputs_cor.mat'),'fval_inv');
cor_fval_ageinv = COR_AGEINV.fval_inv;
% COR_INCOR_TTEST = load(fullfile(maindir,'two_sample_ttest_cor_incor.mat'),'t');
% cor_incor_tval = COR_INCOR_TTEST.t;
% COR_ERRCOR_TTEST = load(fullfile(maindir,'two_sample_ttest_cor_errcor.mat'),'t');
% cor_errcor_tval = COR_ERRCOR_TTEST.t;

% Names of files to be created
cor_permtfcename = 'permutationTFCE_correct.mat';
% cor_incor_permtfcename = 'permutationTFCE_cor_incor.mat';
% cor_errcor_permtfcename = 'permutationTFCE_cor_errcor.mat';

% Number of permutations
nperm = 1000;

% Get number of channels
nchans = size(cor_tval,1);

% preallocating max permuted tfce score vectors for speed
cor_max_perm_tfce_group = zeros(nchans,nperm);
% cor_max_perm_tfce_age = zeros(nchans,nperm);
cor_max_perm_tfce_ageinv = zeros(nchans,nperm);
% cor_incor_max_perm_tfce_group = zeros(nchans,nperm);
% cor_errcor_max_perm_tfce_group = zeros(nchans,nperm);

c = parcluster('local');
c.NumWorkers = N;
p = gcp('nocreate');
sprintf('Number of pools:\n')
disp(p)
parfor n=1:nperm
    % N = str2double(getenv('NUMBER_OF_PROCESSORS'));
    % disp(N)
    % if isnan(N)
    %     setenv('NUMBER_OF_PROCESSORS','1')
    % end
    sprintf('Iteration Number %d',n)
    tic
    % getting total size of t-values matrix
    cor_totalSize = size(cor_tval,1)*size(cor_tval,2)*size(cor_tval,3);
    % cor_incor_totalSize = size(cor_incor_tval,1)*size(cor_incor_tval,2)*size(cor_incor_tval,3);
    % cor_errcor_totalSize = size(cor_errcor_tval,1)*size(cor_errcor_tval,2)*size(cor_errcor_tval,3);

    % reshaping tfce scores to a 1 x totalSize vector
    cor_reshaped_t = reshape(cor_tval,[1,cor_totalSize]);
    % cor_reshaped_f = reshape(cor_fval_age,[1, cor_totalSize]);
    cor_reshaped_finv = reshape(cor_fval_ageinv,[1,cor_totalSize]);
    % cor_incor_reshaped_t = reshape(cor_incor_tval,[1,cor_incor_totalSize]);
    % cor_errcor_reshaped_t = reshape(cor_errcor_tval,[1,cor_errcor_totalSize]);
    
    % Permute t-value matrix
    cor_perm_t = cor_reshaped_t(randperm(length(cor_reshaped_t)));
    % cor_perm_f = cor_reshaped_f(randperm(length(cor_reshaped_f)));
    cor_perm_finv = cor_reshaped_finv(randperm(length(cor_reshaped_finv)));
    % cor_incor_perm_t = cor_incor_reshaped_t(randperm(length(cor_incor_reshaped_t)));
    % cor_errcor_perm_t = cor_errcor_reshaped_t(randperm(length(cor_errcor_reshaped_t)));

    % Reshape permuted t-values back to original size
    cor_perm_t = reshape(cor_perm_t,[size(cor_tval,1) size(cor_tval,2) size(cor_tval,3)]);
    % cor_perm_f = reshape(cor_perm_f,[size(cor_tval,1) size(cor_tval,2) size(cor_tval,3)]);
    cor_perm_finv = reshape(cor_perm_finv, [size(cor_tval,1) size(cor_tval,2) size(cor_tval,3)]);
    % cor_incor_perm_t = reshape(cor_incor_perm_t,[size(cor_incor_tval,1) size(cor_incor_tval,2) size(cor_incor_tval,3)]);
    % cor_errcor_perm_t = reshape(cor_errcor_perm_t,[size(cor_errcor_tval,1) size(cor_errcor_tval,2) size(cor_errcor_tval,3)]);

    % updating waitbar only when running on rhea
%     if strcmp(hostname(1:end-1),'rhea.wpic.upmc.edu')
%         waitbar(n/nperm,f,sprintf('%d',n))
%     end

    % Run TFCE on permuted t-values/ F-values
    % limo_tfce inputs:
    %   type: 3 (for 3D data) where data is channels x freqs x times
    %   data: 3D data
    %   channeighbmatrix: logical matrix indicating neighboring channels
    %   updatebar: 0 or 1 to display TFCE thresholding progress
    %               0 = do not show progress bar
    %               1 = show progress bar
    cor_perm_tfce_score_group = limo_tfce(3,cor_perm_t,channeighbmatrix,0);
    % cor_perm_tfce_score_age = limo_tfce(3,cor_perm_f,channeighbmatrix,0);
    cor_perm_tfce_score_ageinv = limo_tfce(3,cor_perm_finv,channeighbmatrix,0);
    % cor_incor_perm_tfce_score_group = limo_tfce(3,cor_incor_perm_t,channeighbmatrix,0);
    % cor_errcor_perm_tfce_score_group = limo_tfce(3,cor_errcor_perm_t,channeighbmatrix,0);
    
    oneiteration_group = max(cor_perm_tfce_score_group,[],[2 3]);
    oneiteration_ageinv = max(cor_perm_tfce_score_ageinv,[],[2 3]);

    % Save maximum TFCE score over all channels
    cor_max_perm_tfce_group(:,n) = oneiteration_group;
    % cor_max_perm_tfce_age(:,n) = max(cor_perm_tfce_score_age,[],[2,3]);
    cor_max_perm_tfce_ageinv(:,n) = oneiteration_ageinv;
    toc
    savename = fullfile(tfcepath,'permutations_idv',sprintf('%s_%d.mat',cor_permtfcename,n));
    save(savename,'oneiteration_group','oneiteration_ageinv','-fromstruct')
    % cor_incor_max_perm_tfce_group(:,n) = max(cor_incor_perm_tfce_score_group,[],[2 3]);
    % cor_errcor_max_perm_tfce_group(:,n) = max(cor_errcor_perm_tfce_score_group,[],[2 3]);
end
save(fullfile(tfcepath,cor_permtfcename),'cor_max_perm_tfce_ageinv','cor_max_perm_tfce_group')
% save(fullfile(tfcepath,cor_permtfcename),'cor_max_perm_tfce_group','cor_max_perm_tfce_age','cor_max_perm_tfce_ageinv')
% save(fullfile(tfcepath,cor_incor_permtfcename),'cor_incor_max_perm_tfce_group')
% save(fullfile(tfcepath,cor_errcor_permtfcename),'cor_errcor_max_perm_tfce_group')

