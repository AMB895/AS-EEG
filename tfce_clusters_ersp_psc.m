%% Significant Clusters on ERSP
clear; close all
% adding necessary paths
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/tools/')
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/')
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/')
addpath('/ocean/projects/soc230004p/shared/antisaccade_eeg/code/erspFunctions/')
%% Load in correct, vgs, and error data
% Correct AS Trials
% load in ersp data
corerspstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/CorAS/CorASersp.mat');
corersp = corerspstruct.allersp;
times = corerspstruct.preptimes;
freqs = corerspstruct.freqs;
coridmatstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/CorAS/CorASidmatrix.mat');
coridmat = coridmatstruct.idmatrix;
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
vgserspstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/VGS/VGSersp.mat');
vgsersp = vgserspstruct.allersp;
times = vgserspstruct.preptimes;
freqs = vgserspstruct.freqs;
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
errcorerspstruct = load('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ErrCorAS/ErrCorASersp.mat');
errcorersp = errcorerspstruct.allersp;
times = errcorerspstruct.preptimes;
freqs = errcorerspstruct.freqs;
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

% big cell for all ersp data SKIPPING ERROR FOR NOW
% dimensions: 3x2 columns=[frow prow] rows=[cor vgs errcor]
allerspcell = {corersp_viable_frow corersp_viable_prow;vgsersp_frow vgsersp_prow;errcorersp_viable_frow errcorersp_viable_prow};
allidmatcell = {coridmat_viable;vgsidmat;errcoridmat_viable};

% look at models with linear age and inverse age
agetypestrs = ["linage","invage"];

%% looping through linear age and inverse age
for agetype = 1:length(agetypestrs)
%% looping through F row and P rows
for currentRow = 1:length(rownames)
    %% looping through each trial type
    for currentTrialType = 1:length(trialtypenames)
        % define current electrode row and trial type
        currentERSPdata = allerspcell{currentTrialType,currentRow};
        % define current id matrix 
        currentidmat = allidmatcell{currentTrialType};
        % define current age type
        currentagetype = agetypestrs(agetype);
        if strcmp(currentagetype,"linage")
            % set up table for Group Act lme and linear age lme
            % id should be a categorical variable
            T = table('Size',[size(currentidmat,1) 4],'VariableTypes',{'categorical','double','double','double'},...
                'VariableNames',{'id','visit','age','ersp'});
            T.id = currentidmat(:,1);
            T.age = currentidmat(:,3);
            T.visit = currentidmat(:,4);
        elseif strcmp(currentagetype,"invage")
            % set up table for Group Act lme and inverse age lme
            T = table('Size',[size(currentidmat,1) 4],'VariableTypes',{'categorical','double','double','double'},...
                'VariableNames',{'id','visit','age','ersp'});
            T.id = currentidmat(:,1);
            T.age = 1./currentidmat(:,3);
            T.visit = currentidmat(:,4);
        end
        
        %% Group Activation
        % define save name
        % only compute once (no age effects in model)
        groupact_savename = sprintf('%s_groupAct_%s.mat',trialtypenames(currentTrialType),rownames(currentRow));
        
        % check if group activation clusters are already computed
        if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',groupact_savename),'file')
            fprintf('skipping; Computed %s group activation clusters\n',trialtypenames(currentTrialType))
        else
            % updating progress
            fprintf('Group Activation for %s for %s\n',trialtypenames(currentTrialType),rownames(currentRow))
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
        
        %% Age effects (linear or inverse)
        ageeffects_savename = sprintf('%s_%seffect_%s.mat',trialtypenames(currentTrialType),currentagetype,rownames(currentRow));
        
        % check if inverse age clusters are already computed
        if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',ageeffects_savename),'file')
            fprintf('skipping; Computed %s %s effect clusters\n',trialtypenames(currentTrialType),currentagetype)
        else
            % updating progress
            fprintf('Age Effects for %s for %s\n',trialtypenames(currentTrialType),rownames(currentRow))
            % Run linear mixed effects 
            [b,t] = calc_ersp_ageeffects_psc(T,currentERSPdata,numTimes,numFreqs);
            % flip sign of b.age and t.age if age type is inverse
            if strcmp(currentagetype,"invage")
                b.age = -1.*b.age;
                t.age = -1.*t.age;
            end
            % TFCE on t-values
            tfcescores.age = limo_tfce(2,t.age,[]);
            tfcescores.intercept = limo_tfce(2,t.intercept,[]);
            % Permute tfce scores
            permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
            permtfcescores.intercept = calc_perm_tfce_2d(t.intercept,1000);
            % Find significant clusters
            sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age,permtfcescores.age,0.01);
            sigmask.intercept = calc_thres_mask_tfclusters_2d(tfcescores.intercept,permtfcescores.intercept,0.01);
            % save outputs
            savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',ageeffects_savename);
            save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
            % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
            clear t b tfcescores permtfcescores sigmask
        end
    end
    
    %% Models looking at effect of trial type
    % determine age type
    if strcmp(currentagetype,"linage")
        % setting up tables
        Tcor = table('Size',[size(coridmat_viable,1) 5],'VariableTypes',{'categorical','double','double','categorical','double'},...
            'VariableNames',{'id','visit','age','trialtype','ersp'});
        Tcor.id = coridmat_viable(:,1);
        Tcor.visit = coridmat_viable(:,4);
        Tcor.age = coridmat_viable(:,3);
        Tcor.trialtype = repmat("cor",size(coridmat_viable,1),1);

        Tvgs = table('Size',[size(vgsidmat,1) 5],'VariableTypes',{'categorical','double','double','categorical','double'},...
            'VariableNames',{'id','visit','age','trialtype','ersp'});
        Tvgs.id = vgsidmat(:,1);
        Tvgs.visit = vgsidmat(:,4);
        Tvgs.age = vgsidmat(:,3);
        Tvgs.trialtype = repmat("vgs",size(vgsidmat,1),1);

        Terrcor = table('Size',[size(errcoridmat_viable,1) 5],'VariableTypes',{'categorical','double','double','categorical','double'},...
            'VariableNames',{'id','visit','age','trialtype','ersp'});
        Terrcor.id = errcoridmat_viable(:,1);
        Terrcor.visit = errcoridmat_viable(:,4);
        Terrcor.age = errcoridmat_viable(:,3);
        Terrcor.trialtype = repmat("errcor",size(errcoridmat_viable,1),1);
    
    elseif strcmp(currentagetype,"invage")
        % setting up tables
        Tcor = table('Size',[size(coridmat_viable,1) 5],'VariableTypes',{'categorical','double','double','categorical','double'},...
            'VariableNames',{'id','visit','age','trialtype','ersp'});
        Tcor.id = coridmat_viable(:,1);
        Tcor.visit = coridmat_viable(:,4);
        Tcor.age = 1./coridmat_viable(:,3);
        Tcor.trialtype = repmat("cor",size(coridmat_viable,1),1);

        Tvgs = table('Size',[size(vgsidmat,1) 5],'VariableTypes',{'categorical','double','double','categorical','double'},...
            'VariableNames',{'id','visit','age','trialtype','ersp'});
        Tvgs.id = vgsidmat(:,1);
        Tvgs.visit = vgsidmat(:,4);
        Tvgs.age = 1./vgsidmat(:,3);
        Tvgs.trialtype = repmat("vgs",size(vgsidmat,1),1);

        Terrcor = table('Size',[size(errcoridmat_viable,1) 5],'VariableTypes',{'categorical','double','double','categorical','double'},...
            'VariableNames',{'id','visit','age','trialtype','ersp'});
        Terrcor.id = errcoridmat_viable(:,1);
        Terrcor.visit = errcoridmat_viable(:,4);
        Terrcor.age = 1./errcoridmat_viable(:,3);
        Terrcor.trialtype = repmat("errcor",size(errcoridmat_viable,1),1);
    end
    
    %% Effect of trial type on group activation (not controlling for age)
    % Correct AS vs. VGS
    Tcorvgs = [Tcor;Tvgs];
    % checking if file exists
    % only run once b/c no age in model
    corvgs_groupact_savename = sprintf('CorASVGS_trialTypeGroupAct_%s.mat',rownames(currentRow));
    if exist(corvgs_groupact_savename,'file')
        fprintf('skipping; Computed correct AS vs. VGS group activation not controlling for age\n')
    else
        % updating progress
        fprintf('Trial type effect for Correct AS vs. VGS (not controlling for age) for %s\n',rownames(currentRow))
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
    if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_groupact_savename),'file')
        fprintf('skipping; Computed Correct vs. Error Corrected AS group activation not controlling for age\n')
    else
        % updating progress
        fprintf('Trial type effect for Correct vs. Error Corrected AS (not controlling for age) for %s\n',rownames(currentRow))
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
    corvgs_groupactctrl4age_savename = sprintf('CorASVGS_%s_trialTypeGroupAct_%s.mat',currentagetype,rownames(currentRow));
    if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corvgs_groupactctrl4age_savename),'file')
        fprintf('skipping; Computed correct AS vs. VGS group activation controlling for age\n')
    else
        % updating progress
        fprintf('Trial type effect for Correct AS vs. VGS (controlling for %s) for %s\n',currentagetype,rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_trialtype_groupact_psc(Tcorvgs,allerspcell{1,currentRow},allerspcell{2,currentRow},numTimes,numFreqs,1);
        % if inverse age flip sign of b.age and t.age
        if strcmp(currentagetype,"invage")
            b.age = -1.*b.age;
            t.age = -1.*t.age;
        end
        % TFCE on t-values
        tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
        tfcescores.age = limo_tfce(2,t.age,[]);
        tfcescores.intercept = limo_tfce(2,t.intercept,[]);
        % Permute t values
        permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
        permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
        permtfcescores.intercept = calc_perm_tfce_2d(t.intercept,1000);
        % Find significant clusters
        sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype, permtfcescores.trialtype,0.01);
        sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age, permtfcescores.age,0.01);
        sigmask.intercept = calc_thres_mask_tfclusters_2d(tfcescores.intercept, permtfcescores.intercept, 0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corvgs_groupactctrl4age_savename);
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
    % Correct AS vs. Error Correct AS
    corerrcor_groupactctrl4age_savename = sprintf('CorASErrCorAS_%s_trialTypeGroupAct_%s.mat',currentagetype,rownames(currentRow));
    % check if file exists
    if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_groupactctrl4age_savename),'file')
        fprintf('skipping; Computed Correct AS vs. Error Corrected AS group activation controlling for age\n')
    else
        % updating progress
        fprintf('Trial type for Correct vs. Error Corrected AS (controlling for %s) for %s\n',currentagetype,rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_trialtype_groupact_psc(Tcorerrcor,allerspcell{1,currentRow},allerspcell{3,currentRow},numTimes,numFreqs,1);
        % flip sign of b.age and t.age if inverse age
        if strcmp(currentagetype,"invage")
            b.age = -1.*b.age;
            t.age = -1.*t.age;
        end
        % TFCE on t-values
        tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
        tfcescores.age = limo_tfce(2,t.age,[]);
        tfcescores.intercept = limo_tfce(2, t.intercept,[]);
        % Permute t values
        permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
        permtfcescores.age = calc_perm_tfce_2d(t.age, 1000);
        permtfcescores.intercept = calc_perm_tfce_2d(t.intercept, 1000);
        % Find significant clusters
        sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype, permtfcescores.trialtype,0.01);
        sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age, permtfcescores.age, 0.01);
        sigmask.intercept = calc_thres_mask_tfclusters_2d(tfcescores.intercept, permtfcescores.intercept, 0.01);
        % save outputs
        savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_groupactctrl4age_savename);
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
    %% LME for effects of inverse age, trial type and interaction
    % Correct AS vs. VGS
    corvgs_fullLME_savename = sprintf('CorASVGS_fullLME_%s_%s',currentagetype,rownames(currentRow));
    % check if file exists
    if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corvgs_fullLME_savename),'file')
        fprintf('skipping; Computed full LME for Correct AS vs. VGS\n')
    else
        % updating progress
        fprintf('Computing full LME for correct AS vs. VGS for %s\n',rownames(currentRow))
        % Run linear mixed effects model
        [b,t] = calc_ersp_fullLME_psc(Tcorvgs,allerspcell{1,currentRow},allerspcell{2,currentRow},numTimes,numFreqs);
        % if inverse age flip sign of b.age and t.age
        if strcmp(currentagetype,"invage")
            b.age = -1.*b.age;
            t.age = -1.*t.age;
        end
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
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
    
     % Correct AS vs. Error correct AS skipping for now
    corerrcor_fullLME_savename = sprintf('CorASErrCorAS_fullLME_%s',rownames(currentRow));
    % check if file exists
    if exist(sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',corerrcor_fullLME_savename),'file')
        fprintf('skipping; Computed full LME for Correct AS vs. Errcor Corrected AS\n')
    else
       % updating progress
        fprintf('Computing full LME for correct AS vs. error corrected AS for %s\n',rownames(currentRow)) 
        % Run linear mixed effects model
        [b,t] = calc_ersp_fullLME_psc(Tcorerrcor,allerspcell{1,currentRow},allerspcell{3,currentRow},numTimes,numFreqs);
        % if inverse age flip sign of b.age and t.age
        if strcmp(currentagetype,"invage")
            b.age = -1.*b.age;
            t.age = -1.*t.age;
        end
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
        save(savepath,'t','b','tfcescores','permtfcescores','sigmask')
        % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
        clear t b tfcescores permtfcescores sigmask
    end
end
end