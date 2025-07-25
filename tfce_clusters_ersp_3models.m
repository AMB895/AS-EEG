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

% get number of times and frequencies to loop through
numTimes = length(times);
numFreqs = length(freqs);

% set up tables
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

% combine tables
Tcorvgs = [Tcor;Tvgs];
Tcorerrcor = [Tcor;Terrcor];

% big cell for all tables
Tcell = {Tcorvgs; Tcorerrcor};

% Trial types
trialtypenames = ["CorASvsVGS","CorvsErrAS"];

% big cell for all ersp data
corvgs_ersp_cell = {corersp_viable_frow; vgsersp_frow};
corerrcor_ersp_cell = {corersp_viable_frow; errcorersp_viable_frow};
allerspcell = {corvgs_ersp_cell;corerrcor_ersp_cell};

%% looping through each correct AS vs. each trial type
for currentTrialType = 1:length(trialtypenames)
    % defining current comparison and respective ersp data
    currentERSPcell = allerspcell{currentTrialType};
    corersp = currentERSPcell{1};
    ctrlersp = currentERSPcell{2};
    T = Tcell{currentTrialType};
    
    %% Model 1: ERSP ~ 1 + Age + TrialType + Age*TrialType + (1 | lunaid)
    % save name
    m1_savename = sprintf("m1_%s.mat",trialtypenames(currentTrialType));
    % run LMER
    [b,t,AIC] = calc_ersp_fullLME_psc(T,corersp,ctrlersp,numTimes,numFreqs);
    
    % TFCE on t-values
    tfcescores.age = limo_tfce(2,t.age,[]);
    tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
    tfcescores.interaction = limo_tfce(2,t.interaction,[]);
    tfcescores.intercept = limo_tfce(2,t.intercept,[]);
    
    % Permute t-values
    permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
    permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
    permtfcescores.interaction = calc_perm_tfce_2d(t.interaction,1000);
    permtfcescores.intercept = calc_perm_tfce_2d(t.intercept,1000);
    
    % Find significant clusters
    sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age,permtfcescores.age,0.01);
    sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype,permtfcescores.trialtype,0.01);
    sigmask.interaction = calc_thres_mask_tfclusters_2d(tfcescores.interaction,permtfcescores.interaction,0.01);
    sigmask.intercept = calc_thres_mask_tfclusters_2d(tfcescores.intercept,permtfcescores.intercept,0.01);
    
    % save outputs
    savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',m1_savename);
    save(savepath,'t','b','tfcescores','permtfcescores','sigmask','AIC')
    
    % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
    clear t b tfcescores permtfcescores sigmask AIC
    
    %% Model 2 ERSP ~ Age:TrialType + (1 | lunaid)
    % save name
    m2_savename = sprintf("m2_%s.mat",trialtypenames(currentTrialType));
    % run LMER
    [b,t,AIC] = calc_ersp_interaction_psc(T,corersp,ctrlersp,numTimes,numFreqs);
    
    % TFCE on t-values
    tfcescores.interaction = limo_tfce(2,t.interaction,[]);
    tfcescores.intercept = limo_tfce(2,t.intercept,[]);
    
    % Permute t-values
    permtfcescores.interaction = calc_perm_tfce_2d(t.interaction,1000);
    permtfcescores.intercept = calc_perm_tfce_2d(t.intercept,1000);
    
    % Find significant clusters
    sigmask.interaction = calc_thres_mask_tfclusters_2d(tfcescores.interaction,permtfcescores.interaction,0.01);
    sigmask.intercept = calc_thres_mask_tfclusters_2d(tfcescores.intercept,permtfcescores.intercept,0.01);
    
    % save outputs
    savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',m2_savename);
    save(savepath,'t','b','tfcescores','permtfcescores','sigmask','AIC')
    
    % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
    clear t b tfcescores permtfcescores sigmask AIC
    
    %% Model 3: ERSP ~ Age + TrialType + (1 | lunaid)
    % save name
    m3_savename = sprintf("m3_%s.mat",trialtypenames(currentTrialType));
    % run LMER
    [b,t,AIC] = calc_ersp_trialtype_groupact_psc(T,corersp,ctrlersp,numTimes,numFreqs,1);
    
    % TFCE on t-values
    tfcescores.age = limo_tfce(2,t.age,[]);
    tfcescores.trialtype = limo_tfce(2,t.trialtype,[]);
    tfcescores.intercept = limo_tfce(2,t.intercept,[]);
    
    % Permute t-values
    permtfcescores.age = calc_perm_tfce_2d(t.age,1000);
    permtfcescores.trialtype = calc_perm_tfce_2d(t.trialtype,1000);
    permtfcescores.intercept = calc_perm_tfce_2d(t.intercept,1000);
    
    % Find significant clusters
    sigmask.age = calc_thres_mask_tfclusters_2d(tfcescores.age,permtfcescores.age,0.01);
    sigmask.trialtype = calc_thres_mask_tfclusters_2d(tfcescores.trialtype,permtfcescores.trialtype,0.01);
    sigmask.intercept = calc_thres_mask_tfclusters_2d(tfcescores.intercept,permtfcescores.intercept,0.01);
    
    % save outputs
    savepath = sprintf('/ocean/projects/soc230004p/shared/antisaccade_eeg/data/ClusterStats/%s',m3_savename);
    save(savepath,'t','b','tfcescores','permtfcescores','sigmask','AIC')
    
    % clear t, b, tfcescores, permtfcescores, and sigmask for next lme
    clear t b tfcescores permtfcescores sigmask AIC
end
