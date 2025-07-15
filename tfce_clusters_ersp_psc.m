%% Significant Clusters on ERSP
clear; close all
% adding necessary paths
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/tools/')
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/')
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats')
%% Load in correct, vgs, and error data
% Correct AS Trials
% load in ersp data
% corerspstruct = load('/Volumes/Hera/Abby/AS_EEG/ERSPdata/CorASprep/CorASersp.mat');
corerspstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/CorAS/CorASersp.mat');
corersp = corerspstruct.allersp;
times = corerspstruct.preptimes;
freqs = corerspstruct.freqs;
% coridmatstruct = load('/Volumes/Hera/Abby/AS_EEG/ERSPdata/CorASprep/CorASidmatrix.mat');
coridmatstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/CorAS/CorASidmatrix.mat');
coridmat = coridmatstruct.idmatrix;
% load('/Volumes/Hera/Abby/AS_EEG/Results/ErrorLatencyTable_20250320.mat')
load('/ocean/projects/soc230004p/shared/antisaccade_eeg/ErrorLatencyTable_20250320.mat')

% Remove subjects that are not viable
numSubs = size(coridmat,1);
corersp_viable = [];
coridmat_viable = [];
for currentSub = 1:numSubs
    id = coridmat(currentSub,1);
    scandate = coridmat(currentSub,2);
    idx = find(ErrorLatencyTable.LunaID == id & ErrorLatencyTable.ScanDate==scandate);
    subtable = ErrorLatencyTable(idx,:);
    isviable = subtable.Viable==1;
    if isviable
        corersp_viable(end+1,:,:,:) = corersp(currentSub,:,:,:);
        coridmat_viable(end+1,:) = coridmat(currentSub,:,:,:);
    end
end

% Add visit number to corIDmatrix
corIDs = unique(coridmat_viable(:,1));
for currentSub = 1:length(corIDs)
    subIdx = find(coridmat_viable(:,1)==corIDs(currentSub));
    subinfo = coridmat_viable(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       coridmat_viable(subIdx(currentVisit),4) = currentVisit; 
    end
end

% Average ersp across f-row and p-row electrodes
corersp_viable_frow = squeeze(mean(corersp_viable(:,[4 5 6 7 37 38 39 40 41],:,:),2));
corersp_viable_prow = squeeze(mean(corersp_viable(:,[19 20 21 22 23 55 56 57 58 59 31],:,:),2));

% VGS Trials
% load in VGS data
% vgserspstruct = load('/Volumes/Hera/Abby/AS_EEG/ERSPdata/VGSprep/VGSersp.mat');
vgserspstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/VGS/VGSersp.mat');
vgsersp = vgserspstruct.allersp;
times = vgserspstruct.preptimes;
freqs = vgserspstruct.freqs;
% vgsidmatstruct = load('/Volumes/Hera/Abby/AS_EEG/ERSPdata/VGSprep/VGSidmatrix.mat');
vgsidmatstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/VGS/VGSidmatrix.mat');
vgsidmat = vgsidmatstruct.idmatrix;

% add visit number to vgsidmat
vgsIDs = unique(vgsidmat(:,1));
for currentSub = 1:length(vgsIDs)
    subIdx = find(vgsidmat(:,1)==vgsIDs(currentSub));
    subinfo = vgsidmat(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       vgsidmat(subIdx(currentVisit),4) = currentVisit; 
    end
end

% average across f-row
vgsersp_frow = squeeze(mean(vgsersp(:,[4 5 6 7 37 38 39 40 41],:,:),2));
vgsersp_prow = squeeze(mean(vgsersp(:,[19 20 21 22 23 55 56 57 58 59 31],:,:),2));

% Error corrected AS trials
% load in error data
% errcorerspstruct = load('/Volumes/Hera/Abby/AS_EEG/ERSPdata/ErrCorASprep/ErrCorASersp.mat');
errcorerspstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ErrCorAS/ErrCorASersp.mat');
errcorersp = errcorerspstruct.allersp;
times = errcorerspstruct.preptimes;
freqs = errcorerspstruct.freqs;
% errcoridmatstruct = load('/Volumes/Hera/Abby/AS_EEG/ERSPdata/ErrCorASprep/ErrCorASidmatrix.mat');
errcoridmatstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ErrCorAS/ErrCorASidmatrix.mat');
errcoridmat = errcoridmatstruct.idmatrix;
% Remove subjects that are not viable
numSubs = size(errcoridmat,1);
errcorersp_viable = [];
errcoridmat_viable = [];
for currentSub = 1:numSubs
    id = errcoridmat(currentSub,1);
    scandate = errcoridmat(currentSub,2);
    idx = find(ErrorLatencyTable.LunaID == id & ErrorLatencyTable.ScanDate==scandate);
    subtable = ErrorLatencyTable(idx,:);
    isviable = subtable.Viable==1;
    if isviable
        errcorersp_viable(end+1,:,:,:) = errcorersp(currentSub,:,:,:);
        errcoridmat_viable(end+1,:) = errcoridmat(currentSub,:,:,:);
    end
end

% Add visit number to errcorIDmatrix
errcorIDs = unique(errcoridmat_viable(:,1));
for currentSub = 1:length(errcorIDs)
    subIdx = find(errcoridmat_viable(:,1)==errcorIDs(currentSub));
    subinfo = errcoridmat_viable(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       errcoridmat_viable(subIdx(currentVisit),4) = currentVisit; 
    end
end
% Average ersp across f-row and p-row electrodes
errcorersp_viable_frow = squeeze(mean(errcorersp_viable(:,[4 5 6 7 37 38 39 40 41],:,:),2));
errcorersp_viable_prow = squeeze(mean(errcorersp_viable(:,[19 20 21 22 23 55 56 57 58 59 31],:,:),2));

% get number of times and frequencies to loop through
numTimes = length(times);
numFreqs = length(freqs);

% Row names and ERSP data
rownames = ["Frow","Prow"];

% Trial types
trialtypenames = ["CorAS","VGS","ErrCorAS"];

% big cell for all ersp data
% dimensions: 3x2 columns=[frow prow] rows=[cor vgs errcor]
allerspcell = {corersp_viable_frow corersp_viable_prow;vgsersp_frow vgsersp_prow;errcorersp_viable_frow errcorersp_viable_prow};
allidmatcell = {coridmat_viable;vgsidmat;errcoridmat_viable};

%% First looping through F row and P rows
for currentRow = 1:length(rownames)
    % looping through each trial type
    for currentTrialType = 1:length(trialtypenames)
        % define current electrode row and trial type
        currentERSPdata = allerspcell{currentTrialType,currentRow};
        % define current id matrix 
        currentidmat = allidmatcell{currentTrialType};
        % set up table for Group Act lme and inverse age lme
        T = table('Size',[size(currentidmat,1) 4],'VariableTypes',{'double','double','double','double'},...
            'VariableNames',{'id','visit','invage','ersp'});
        T.id = currentidmat(:,1);
        T.invage = 1./currentidmat(:,3);
        T.visit = currentidmat(:,4);
        
        %% Group Activation
        % define save names
        groupact_savename = sprintf('%s_groupAct_%s.mat',trialtypenames(currentTrialType),rownames(currentRow));
        
        % check if group activation clusters are already computed
        if exist(groupact_savename,'file')
            fprintf('Computed %s group activation clusters\n',trialtypenames(currentTrialType))
        else
            % updating progress
            fprintf('Computing group Activation for %s for %s\n',trialtypenames(currentTrialType),rownames(currentRow))
            % Run linear mixed effects model for group activation
            [b,t] = calc_ersp_groupact_psc(T,currentERSPdata,numTimes,numFreqs);
            % TFCE on t-values
            tfcescores = limo_tfce(2,t,[]);
            % Permute TFCE scores
            permtfcescores = calc_perm_tfce_2d(t,1000);
            % Find significant clusters
            sigmask = calc_thres_mask_tfclusters_2d(tfcescores,permtfcescores,0.01);
            % save outputs
            savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',groupact_savename);
            save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
            % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
            clear t b tfcescores permtfcescores sigmask
        end
        
        %% Inverse Age effects
        inverseage_savename = sprintf('%s_inverseage_%s.mat',trialtypenames(currentTrialType),rownames(currentRow));
        
        % check if inverse age clusters are already computed
        if exist(inverseage_savename,'file')
            fprintf('Computed %s inverse age clusters\n',trialtypenames(currentTrialType))
        else
            % updating progress
            fprintf('Inverse Age Effects for %s for %s\n',trialtypenames(currentTrialType),rownames(currentRow))
            % Run linear mixed effects 
            [b,t] = calc_ersp_invageeffects_psc(T,currentERSPdata,numTimes,numFreqs);
            % TFCE on t-values
            tfcescores.age = limo_tfce(2,t.age,[]);
            tfcescores.groupact_ctrl4age = limo_tfce(2,t.intercept,[]);
            % Permute tfce scores
            permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
            permtfcescores.groupact_ctrl4age = calc_perm_tfce_2d(t.intercept,1000);
            % Find significant clusters
            sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age,permtfcescores.age,0.01);
            sigmask.groupact_ctrl4age = calc_thres_mask_tfclusters_2d(tfcescores.groupact_ctrl4age,permtfcescores.groupact_ctrl4age,0.01);
            % save outputs
            savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',inverseage_savename);
            save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
            % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
            clear t b tfcescores permtfcescores sigmask
        end
    end
    
    %% Models looking at effect of trial type
    % setting up tables
    Tcor = table('Size',[size(coridmat_viable,1) 5],'VariableTypes',{'double','double','double','categorical','double'},...
        'VariableNames',{'id','visit','invage','trialtype','ersp'});
    Tcor.id = coridmat_viable(:,1);
    Tcor.visit = coridmat_viable(:,4);
    Tcor.invage = 1./coridmat_viable(:,3);
    Tcor.trialtype = repmat("cor",size(coridmat_viable,1),1);
    
    Tvgs = table('Size',[size(vgsidmat,1) 5],'VariableTypes',{'double','double','double','categorical','double'},...
        'VariableNames',{'id','visit','invage','trialtype','ersp'});
    Tvgs.id = vgsidmat(:,1);
    Tvgs.visit = vgsidmat(:,4);
    Tvgs.invage = 1./vgsidmat(:,4);
    Tvgs.trialtype = repmat("vgs",size(vgsidmat,1),1);
    
    Terrcor = table('Size',[size(errcoridmat_viable,1) 5],'VariableTypes',{'double','double','double','categorical','double'},...
        'VariableNames',{'id','visit','invage','trialtype','ersp'});
    Terrcor.id = errcoridmat_viable(:,1);
    Terrcor.visit = errcoridmat_viable(:,4);
    Terrcor.invage = 1./errcoridmat_viable(:,3);
    Terrcor.trialtype = repmat("errcor",size(errcoridmat_viable,1),1);
    
    %% Effect of trial type on group activation (not controlling for age)
    % Correct AS vs. VGS
    Tcorvgs = [Tcor;Tvgs];
    % checking if file exists
    corvgs_groupact_savename = sprintf('CorASVGS_trialTypeGroupAct_%s.mat',rownames(currentRow));
    if exist(corvgs_groupact_savename,'file')
        fprintf('Computed correct AS vs. VGS group activation not controlling for age\n')
    else
        % updating progress
        fprintf('Computing effect of trial type for Correct AS vs. VGS (not controlling for age) for %s\n',rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_trialtype_groupact_psc(Tcorvgs,allerspcell{1,currentRow},allerspcell{2,currentRow},numTimes,numFreqs,0);
        % TFCE on t-values
        tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
        % Permute t-values
        permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
        % Find significant clusters
        sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype,permtfcescores.trialtype,0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corvgs_groupact_savename);
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
    % Correct AS vs. Error Correct AS
    Tcorerrcor = [Tcor;Terrcor];
    corerrcor_groupact_savename = sprintf('CorASErrCorAS_trialTypeGroupAct_%s.mat',rownames(currentRow));
    % check if file exists
    if exist(corerrcor_groupact_savename,'file')
        fprintf('Computed group activation by trial type for Correct AS vs. Error Corrected AS not controlling for age\n')
    else
        % updating progress
        fprintf('Computing effect of trial type for Correct AS vs. Error Corrected AS (not controlling for age) for %s\n',rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_trialtype_groupact_psc(Tcorerrcor,allerspcell{1,currentRow},allerspcell{3,currentRow},numTimes,numFreqs,0);
        % TFCE on t-values
        tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
        % Permute t-values
        permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
        % Find significant clusters
        sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype,permtfcescores.trialtype,0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_groupact_savename);
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
    %% Effect of trial type on group activation (controlling for age)
    % Correct AS vs. VGS
    % checking if file exists
    corvgs_groupactctrl4age_savename = sprintf('CorASVGS_trialTypeGroupAct_ctrl4age_%s.mat',rownames(currentRow));
    if exist(corvgs_groupactctrl4age_savename,'file')
        fprintf('Computed correct AS vs. VGS group activation controlling for age\n')
    else
        % updating progress
        fprintf('Computing effect of trial type for Correct AS vs. VGS (controlling for age) for %s\n',rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_trialtype_groupact_psc(Tcorvgs,allerspcell{1,currentRow},allerspcell{2,currentRow},numTimes,numFreqs,1);
        % TFCE on t-values
        tfcescores.trialtype_ctrl4age = limo_tfce(2,t.trialtype,[]);
        % Permute t values
        permtfcescores.trialtype_ctrl4age = calc_perm_tfce_2d(t.trialtype,1000);
        % Find significant clusters
        sigmask.trialtype_ctrl4age = calc_thres_mask_tfclusters_2d(tfcescores.trialtype_ctrl4age,permtfcescores.trialtype_ctrl4age,0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corvgs_groupactctrl4age_savename);
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
    % Correct AS vs. Error Correct AS
    corerrcor_groupactctrl4age_savename = sprintf('CorASErrCorAS_trialTypeGroupAct_ctrl4age_%s.mat',rownames(currentRow));
    % check if file exists
    if exist(corerrcor_groupactctrl4age_savename,'file')
        fprintf('Computed group activation by trial type for Correct AS vs. Error Corrected AS controlling for age\n')
    else
        % updating progress
        fprintf('Computing effect of trial type for Correct AS vs. Error Corrected AS (controlling for age) for %s\n',rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_trialtype_groupact_psc(Tcorerrcor,allerspcell{1,currentRow},allerspcell{3,currentRow},numTimes,numFreqs,1);
        % TFCE on t-values
        tfcescores.trialtype_ctrl4age = limo_tfce(2,t.trialtype,[]);
        % Permute t values
        permtfcescores.trialtype_ctrl4age = calc_perm_tfce_2d(t.trialtype,1000);
        % Find significant clusters
        sigmask.trialtype_ctrl4age = calc_thres_mask_tfclusters_2d(tfcescores.trialtype_ctrl4age,permtfcescores.trialtype_ctrl4age,0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_groupactctrl4age_savename);
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
    %% LME for effects of inverse age, trial type and interaction
    % Correct AS vs. VGS
    corvgs_fullLME_savename = sprintf('CorASVGS_fullGLM_%s',rownames(currentRow));
    % check if file exists
    if exist(corvgs_fullLME_savename,'file')
        fprintf('Computed full LME for Correct AS vs. VGS\n')
    else
        % updating progress
        fprintf('Computing full LME for correct AS vs. VGS for %s\n',rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_fullLME_psc(Tcorvgs,allerspcell{1,currentRow},allerspcell{2,currentRow},numTimes,numFreqs);
        % TFCE on t-values
        tfcescores.age = limo_tfce(2,t.age,[]);
        tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
        tfcescores.interaction = limo_tfce(2,t.interaction,[]);
        % Permute tfce
        permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
        permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
        permtfcescores.interaction = calc_perm_tfce_2d(t.interaction,1000);
        % Find significant clusters
        sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age,permtfcescores.age,0.01);
        sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype,permtfcescores.trialtype,0.01);
        sigmask.interaction = calc_thres_mask_tfclusters_2d(tfcescores.interaction,permtfcescores.interaction,0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corvgs_fullLME_savename);
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
     % Correct AS vs. Error correct AS
    corerrcor_fullLME_savename = sprintf('CorASErrCorAS_fullGLM_%s',rownames(currentRow));
    % check if file exists
    if exist(corerrcor_fullLME_savename,'file')
        fprintf('Computed full LME for Correct AS vs. Errcor Corrected AS\n')
    else
       % updating progress
        fprintf('Computing full LME for correct AS vs. error corrected AS for %s\n',rownames(currentRow)) 
        % Run linear mixed effects model
        [b,t] = calc_ersp_fullLME_psc(Tcorerrcor,allerspcell{1,currentRow},allerspcell{3,currentRow},numTimes,numFreqs);
        % TFCE on t-values
        tfcescores.age = limo_tfce(2,t.age,[]);
        tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
        tfcescores.interaction = limo_tfce(2,t.interaction,[]);
        % Permute tfce
        permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
        permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
        permtfcescores.interaction = calc_perm_tfce_2d(t.interaction,1000);
        % Find significant clusters
        sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age,permtfcescores.age,0.01);
        sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype,permtfcescores.trialtype,0.01);
        sigmask.interaction = calc_thres_mask_tfclusters_2d(tfcescores.interaction,permtfcescores.interaction,0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_fullLME_savename);
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
end
