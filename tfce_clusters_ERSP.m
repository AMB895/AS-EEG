%% Significant Clusters on ERSP data
% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')
addpath('/Volumes/Hera/Abby/AS_EEG/erspFunctions/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/')

%% Correct AS trials
clear; close all;
% load in ERSP data
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable_20250320.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/corERSPdata.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/corIDmatrix.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/cormissingdata.mat')

%% First remove participants that are not viable
numSubs = size(corIDmatrix,1);

corerspdata_viable=[];
corIDmatrix_viable=[];
for currentSub = 1:numSubs
    id = corIDmatrix(currentSub,1);
    scandate = corIDmatrix(currentSub,2);
    idx = find(ErrorLatencyTable.LunaID == id & ErrorLatencyTable.ScanDate==scandate);
    subtable = ErrorLatencyTable(idx,:);
    isviable = subtable.Viable==1;
    if isviable
        corerspdata_viable(end+1,:,:,:) = corerspdata(currentSub,:,:,:);
        corIDmatrix_viable(end+1,:) = corIDmatrix(currentSub,:);
    end
end

% Add visit number to corIDmatrix
corIDs = unique(corIDmatrix_viable(:,1));
for currentSub = 1:length(corIDs)
    subIdx = find(corIDmatrix_viable(:,1)==corIDs(currentSub));
    subinfo = corIDmatrix_viable(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       corIDmatrix_viable(subIdx(currentVisit),4) = currentVisit; 
    end
end

% Average across F-row
corerspdata_viable_frow = squeeze(mean(corerspdata_viable(:,[4 5 6 7 37 38 39 40],:,:),2));

%% Group Activation- Correct Trials
numTimes = length(times);
numFreqs = length(freqs);

T = table('Size',[size(corIDmatrix_viable,1) 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'id','visit','power'});
T.id = corIDmatrix_viable(:,1);
T.visit = corIDmatrix_viable(:,4);

if exist('corGroupActClusters.mat','file')
    fprintf('Computed correct group activation clusters; loading\n')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corGroupActClusters.mat')
else
    % Run group activation linear model
    [b_cor_groupact,t_cor_groupact,p_cor_groupact] = calc_ersp_groupact(T,corerspdata_viable_frow,numTimes,numFreqs);
    % TFCE on t-values
    tfcescores_cor_groupact = limo_tfce(2,t_cor_groupact,[]);
    % Permute TFCE scores
    permtfcescores_cor_groupact = calc_perm_tfce_2d(t_cor_groupact,1000);
    % Find significant clusters
    [~,mask_cor_groupact] = calc_thres_mask_tfclusters_2d(tfcescores_cor_groupact,permtfcescores_cor_groupact,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corGroupActClusters.mat','t_cor_groupact','b_cor_groupact','p_cor_groupact','tfcescores_cor_groupact','mask_cor_groupact')
end

% Plot t-values with mask
figure;
surf(times,freqs,mask_cor_groupact.*t_cor_groupact,'EdgeColor','none')
C = colorbar;
C.Label.String = 't stat';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Correct Trials Group Activation')
subtitle('F-row electrodes')

%% Inverse Age Effects
T = table('Size',[size(corIDmatrix_viable,1) 4],'VariableTypes',{'double','double','double','double'},...
    'VariableNames',{'id','visit','invage','power'});
T.id = corIDmatrix_viable(:,1);
T.visit = corIDmatrix_viable(:,4);
T.invage = 1./corIDmatrix_viable(:,3);

if exist('corAgeEffectClusters.mat','file')
    fprintf('Computed correct inverse age effect clusters; loading\n')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corAgeEffectClusters.mat')
else
    % Run linear regression with inverse age
    [b_cor_invage,t_cor_invage,p_cor_invage]=calc_ersp_invageeffects(T,corerspdata_viable_frow,numTimes,numFreqs);
    % TFCE on F-values
    tfcescores_cor_invage = limo_tfce(2,t_cor_invage,[]);
    % Permute TFCE scores
    permtfcescores_cor_invage = calc_perm_tfce_2d(t_cor_invage,1000);
    % Find significant clusters
    [~,mask_cor_invage]  = calc_thres_mask_tfclusters_2d(tfcescores_cor_invage,permtfcescores_cor_invage,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corAgeEffectClusters.mat','t_cor_invage','b_cor_invage','p_cor_invage','tfcescores_cor_invage','mask_cor_invage')
end

% Plot age coefficients with mask
figure;
surf(times,freqs,mask_cor_invage.*b_cor_invage,'EdgeColor','none')
C = colorbar;
C.Label.String = 'Age Coefficient (slope)';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Correct Trials Age Effects')
subtitle('F-row electrodes')
%% Error Corrected AS trials
% load in ERSP data
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/errcorERSPdata.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/errcorIDmatrix.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/errcormissingdata.mat')

%% First remove participants that are not viable
numSubs = size(errcorIDmatrix,1);

errcorerspdata_viable=[];
errcorIDmatrix_viable=[];
for currentSub = 1:numSubs
    id = errcorIDmatrix(currentSub,1);
    scandate = errcorIDmatrix(currentSub,2);
    idx = find(ErrorLatencyTable.LunaID == id & ErrorLatencyTable.ScanDate==scandate);
    subtable = ErrorLatencyTable(idx,:);
    isviable = subtable.Viable==1;
    if isviable
        errcorerspdata_viable(end+1,:,:,:) = errcorerspdata(currentSub,:,:,:);
        errcorIDmatrix_viable(end+1,:) = errcorIDmatrix(currentSub,:);
    end
end

% Add visit number to errcorIDmatrix
errcorIDs = unique(errcorIDmatrix_viable(:,1));
for currentSub = 1:length(corIDs)
    subIdx = find(errcorIDmatrix_viable(:,1)==corIDs(currentSub));
    subinfo = errcorIDmatrix_viable(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       errcorIDmatrix_viable(subIdx(currentVisit),4) = currentVisit; 
    end
end
% Average across F-row
errcorerspdata_viable_frow = squeeze(mean(errcorerspdata_viable(:,[4 5 6 7 37 38 39 40],:,:),2));

%% Group Activation- Error Correct Trials
numTimes = length(times);
numFreqs = length(freqs);

T = table('Size',[size(errcorIDmatrix_viable,1) 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'id','visit','power'});
T.id = errcorIDmatrix_viable(:,1);
T.visit = errcorIDmatrix_viable(:,4);
if exist('errcorGroupActClusters.mat','file')
    fprintf('Computed error corrected group activation clusters; loading\n')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/errcorGroupActClusters.mat')
else
    % Run group activation linear model
    [b_errcor_groupact,t_errcor_groupact,p_errcor_groupact] = calc_ersp_groupact(T,errcorerspdata_viable_frow,numTimes,numFreqs);
    % TFCE on t-values
    tfcescores_errcor_groupact = limo_tfce(2,t_errcor_groupact,[]);
    % Permute TFCE scores
    permtfcescores_errcor_groupact = calc_perm_tfce_2d(t_errcor_groupact,1000);
    % Find significant clusters
    [~,mask_errcor_groupact] = calc_thres_mask_tfclusters_2d(tfcescores_errcor_groupact,permtfcescores_errcor_groupact,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/errcorGroupActClusters.mat','t_errcor_groupact','b_errcor_groupact','p_errcor_groupact','tfcescores_errcor_groupact','mask_errcor_groupact')
end

% Plot t-values with mask
figure;
surf(times,freqs,mask_errcor_groupact.*t_errcor_groupact,'EdgeColor','none')
C = colorbar;
C.Label.String = 't stat';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Error Trials Group Activation')
subtitle('F-row electrodes')

%% Inverse Age Effects- Error Correct Trials
T = table('Size',[size(errcorIDmatrix_viable,1) 4],'VariableTypes',{'double','double','double','double'},...
    'VariableNames',{'id','visit','invage','power'});
T.id = errcorIDmatrix_viable(:,1);
T.visit = errcorIDmatrix_viable(:,4);
T.invage = 1./errcorIDmatrix_viable(:,3);

if exist('errcorAgeEffectClusters.mat','file')
    fprintf('Computed error corrected age effect clusters; loading\n')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/errcorAgeEffectClusters.mat')
else
    % Run linear regression with inverse age
    [b_errcor_invage,t_errcor_invage,p_errcor_invage]=calc_ersp_invageeffects(T,errcorerspdata_viable_frow,numTimes,numFreqs);
    % TFCE on F-values
    tfcescores_errcor_invage = limo_tfce(2,t_errcor_invage,[]);
    % Permute TFCE scores
    permtfcescores_errcor_invage = calc_perm_tfce_2d(t_errcor_invage,1000);
    % Find significant clusters
    [~,mask_errcor_invage]  = calc_thres_mask_tfclusters_2d(tfcescores_errcor_invage,permtfcescores_errcor_invage,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/errcorAgeEffectClusters.mat','t_errcor_invage','b_errcor_invage','p_errcor_invage','tfcescores_errcor_invage','mask_errcor_invage')
end
% Plot age coefficients with mask
figure;
surf(times,freqs,mask_errcor_invage.*b_errcor_invage,'EdgeColor','none')
C = colorbar;
C.Label.String = 'Age Coefficient (slope)';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Error Trials Age Effects')
subtitle('F-row electrodes')
%% VGS Trials
% load in ERSP data
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs/vgsERSPdata.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs/vgsIDmatrix.mat')
load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs/vgsmissingdata.mat')

% Add visit number to vgsIDmatrix
vgsIDs = unique(vgsIDmatrix(:,1));
for currentSub = 1:length(vgsIDs)
    subIdx = find(vgsIDmatrix(:,1)==vgsIDs(currentSub));
    subinfo = vgsIDmatrix(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       vgsIDmatrix(subIdx(currentVisit),4) = currentVisit; 
    end
end
% Average across F-row
vgserspdata_frow = squeeze(mean(vgserspdata(:,[4 5 6 7 37 38 39 40],:,:),2));
%% Group Activation- VGS Trials
numTimes = length(times);
numFreqs = length(freqs);

T = table('Size',[size(vgsIDmatrix,1) 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'id','visit','power'});
T.id = vgsIDmatrix(:,1);
T.visit = vgsIDmatrix(:,4);

if exist('vgsGroupActClusters.mat','file')
    fprintf('Computed vgs group activation clusters; loading\n')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgsGroupActClusters.mat')
else
    % Run group activation linear model
    [b_vgs_groupact,t_vgs_groupact,p_vgs_groupact] = calc_ersp_groupact(T,vgserspdata_frow,numTimes,numFreqs);
    % TFCE on t-values
    tfcescores_vgs_groupact = limo_tfce(2,t_vgs_groupact,[]);
    % Permute TFCE scores
    permtfcescores_vgs_groupact = calc_perm_tfce_2d(t_vgs_groupact,1000);
    % Find significant clusters
    [~,mask_vgs_groupact] = calc_thres_mask_tfclusters_2d(tfcescores_vgs_groupact,permtfcescores_vgs_groupact,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgsGroupActClusters.mat','t_vgs_groupact','b_vgs_groupact','p_vgs_groupact','tfcescores_vgs_groupact','mask_vgs_groupact')
end

% Plot t-values with mask
figure;
surf(times,freqs,mask_vgs_groupact.*t_vgs_groupact,'EdgeColor','none')
C = colorbar;
C.Label.String = 't stat';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('VGS Trials Group Activation')
subtitle('F-row electrodes')
%% Inverse Age Effects- VGS Trials
T = table('Size',[size(vgsIDmatrix,1) 4],'VariableTypes',{'double','double','double','double'},...
    'VariableNames',{'id','visit','invage','power'});
T.id = vgsIDmatrix(:,1);
T.visit = vgsIDmatrix(:,4);
T.invage = 1./vgsIDmatrix(:,3);

if exist('vgsAgeEffectClusters.mat','file')
    fprintf('Computed vgs age effect clusters; loading')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgsAgeEffectClusters.mat')
else
% Run linear regression with inverse age
[b_vgs_invage,t_vgs_invage,p_vgs_invage]=calc_ersp_invageeffects(T,vgserspdata_frow,numTimes,numFreqs);
% TFCE on F-values
tfcescores_vgs_invage = limo_tfce(2,t_vgs_invage,[]);
% Permute TFCE scores
permtfcescores_vgs_invage = calc_perm_tfce_2d(t_vgs_invage,1000);
% Find significant clusters
[~,mask_vgs_invage]  = calc_thres_mask_tfclusters_2d(tfcescores_vgs_invage,permtfcescores_vgs_invage,0.05);
% save oututs
save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgsAgeEffectClusters.mat','t_vgs_invage','b_vgs_invage','p_vgs_invage','tfcescores_vgs_invage','mask_vgs_invage')
end

% Plot f-values with mask
figure;
surf(times,freqs,mask_vgs_invage.*b_vgs_invage,'EdgeColor','none')
C = colorbar;
C.Label.String = 'Age Coefficient (slope)';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('VGS Trials Age Effects')
subtitle('F-row electrodes')

%% Correct vs. VGS Trials
% Match IDs from correct AS and VGS trials
corvgsIDmatrix = [];
vgserspdata_corvgs = [];
corerspdata_corvgs = [];
for idxVGS = 1:size(vgsIDmat,1)
   currentSub = vgsIDmat(idxVGS,1);
   currentDate = vgsIDmat(idxVGS,2);
   idxAS = find(corIDmat(:,1)==currentSub & corIDmat(:,2) == currentDate);
   if ~isempty(idxAS) % there is a matching ID and date
       corvgsIDmatrix(end+1,:) = vgsIDmat(idxVGS,:);
       vgserspdata_corvgs(end+1,:,:) = vgserspdata_frow(idxVGS,:,:);
       corerspdata_corvgs(end+1,:,:) = corerspdata_viable_frow(idxAS,:,:);
   end
end

T = table('Size',[size(corvgsIDmatrix,1) 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'id','visit','power'});
T.id = corvgsIDmatrix(:,1);
T.visit = corvgsIDmatrix(:,4);

if exist('corvgsGroupActClusters.mat','file')
    fprintf('Computed correct vs. vgs group activation; loading')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corvgsGroupActClusters.mat')
else
    % Difference ERSP data to run linear model on for activation
    differspdata_corvgs = corerspdata_corvgs - vgserspdata_corvgs;
    % Run group activation linear model
    [b_corvgs_groupact,t_corvgs_groupact,p_corvgs_groupact] = calc_ersp_groupact(T,differspdata_corvgs,numTimes,numFreqs);
    % TFCE on t-values
    tfcescores_corvgs_groupact = limo_tfce(2,t_corvgs_groupact,[]);
    % Permute TFCE scores
    permtfcescores_corvgs_groupact = limo_tfce(2,t_corvgs_groupact,[]);
    % find significant clusters
    [~,mask_corvgs_groupact] = calc_thres_mask_tfclusters_2d(tfcescores_corvgs_groupact,permtfcescores_corvgs_groupact,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corvgsGroupActClusters.mat','b_corvgs_groupact','t_corvgs_groupact','_p_corvgs_groupact','tfcescores_corvgs_groupact','mask_corvgs_groupact')
end

% Plot t-values with mask
figure;
surf(times,freqs,mask_corvgs_groupact.*t_corvgs_groupact,'EdgeColor','none')
C = colorbar;
C.Label.String = 't stat';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Correct AS vs. VGS Trials Group Activation')
subtitle('F-row electrodes')

%% Inverse Age effects, Trial Type effects and Interaction for Correct vs. VGS
% Power ~ TrialType + InvAge + TrialType*InvAge + (1 | ID) @ each time-frequency point
% set up table for each time-frequency point
T = table('Size',[2*size(corvgsIDmatrix,1) 5],'VariableTypes',{'double','double','double','categorical','double'},...
    'VariableNames',{'id','visit','invage','trialtype','power'});
T.id = [corvgsIDmatrix(:,1);corvgsIDmatrix(:,1)];
T.visit = [corvgsIDmatrix(:,4);corvgsIDmatrix(:,4)];
T.invage = [1./corvgsIDmatrix(:,3); 1./corvgsIDmatrix(:,3)];
T.trialtype =[repmat("cor",size(corvgsIDmatrix,1),1);repmat("vgs",size(corvgsIDmatrix,1),1)];
T.trialtype = categorical(T.trialtype);

% Run linear mixed model
[b_corvgs_mixedmodel,t_corvgs_mixedmodel,p_corvgs_mixedmodel] = calc_ersp_linearmixedmodel(T,corerspdata_corvgs,vgserspdata_corvgs,numTimes,numFreqs);

%% Correct vs. Error Trials
% Match IDs from correct AS and erroc corrected AS trials
corerrcorIDmat = [];
corerspdata_corerrcor = [];
errcorerspdata_corerrcor = [];
for idxCOR = 1:size(corIDmatrix_viable,1)
   currentSub = corIDmatrix_viable(idxCOR,1);
   currentDate = corIDmatrix_viable(idxCOR,2);
   idxERRCOR = find(errcorIDmat(:,1)==currentSub & errcorIDmat(:,2) == currentDate);
   if ~isempty(idxERRCOR) % there is a matching ID and date
       corerrcorIDmat(end+1,:) = corIDmatrix_viable(idxCOR,:);
       corerspdata_corerrcor(end+1,:,:) = corerspdata_viable_frow(idxCOR,:,:);
       errcorerspdata_corerrcor(end+1,:,:) = errcorerspdata_viable_frow(idxERRCOR,:,:);
   end
end
T = table('Size',[size(corerrcorIDmat,1) 3],'VariableTypes',{'double','double','double'},...
    'VariableNames',{'id','visit','power'});
T.id = corerrcorIDmat(:,1);
T.visit = corerrcorIDmat(:,4);

if exist('corerrcorGroupActClusters.mat','file')
    fprintf('Computed correct vs. error corrected group activation clusters; loading')
    load('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corerrcorGroupActClusters.mat')
else
    % Difference ERSP data to run linear model on for activation
    differspdata_corerrcor = corerspdata_corerrcor - errcorerspdata_corerrcor;
    % Run group activation linear model
    [b_corerrcor_groupact,t_corerrcor_groupact,~] = calc_ersp_groupact(T,differspdata_corerrcor,numTimes,numFreqs);
    % TFCE on t-values
    tfcescores_corerrcor_groupact = limo_tfce(2,t_corerrcor_groupact,[]);
    % Permute TFCE scores
    permtfcescores_corerrcor_groupact = limo_tfce(2,t_corerrcor_groupact,[]);
    % find significant clusters
    [~,mask_corerrcor_groupact] = calc_thres_mask_tfclusters_2d(tfcescores_corerrcor_groupact,permtfcescores_corerrcor_groupact,0.05);
    % save outputs
    save('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/corerrcorGroupActClusters.mat','b_corerrcor_groupact','t_corerrcor_groupact','p_corerrcor_groupact','tfcescores_corerrcor_groupact','mask_corerrcor_groupact')
end

% Plot t-values with mask
figure;
surf(times,freqs,mask_corerrcor_groupact.*t_corerrcor_groupact,'EdgeColor','none')
C = colorbar;
C.Label.String = 't stat';
view(2)
hold on
xline(0,'--r','LineWidth',1.5);xline(500,'--k','LineWidth',1.5)
xlabel('Time (ms)')
ylabel('Frequency (Hz)')
title('Error AS vs. VGS Trials Group Activation')
subtitle('F-row electrodes')

%% Inverse Age effects, Trial Type effects and Interaction for Correct vs. Error Trials
% Power ~ TrialType + InvAge + TrialType*InvAge + (1 | ID) @ each time-frequency point
% set up table for each time-frequency point
T = table('Size',[2*size(corerrcorIDmat,1) 5],'VariableTypes',{'double','double','double','categorical','double'},...
    'VariableNames',{'id','visit','invage','trialtype','power'});
T.id = [corerrcorIDmat(:,1);corerrcorIDmat(:,1)];
T.visit = [corerrcorIDmat(:,4);corerrcorIDmat(:,4)];
T.invage = [1./corerrcorIDmat(:,3); 1./corerrcorIDmat(:,3)];
T.trialtype =[repmat("cor",size(corerrcorIDmat,1),1);repmat("errcor",size(corerrcorIDmat,1),1)];
T.trialtype = categorical(T.trialtype);

% Run linear mixed model
[b_corerrcor_mixedmodel,t_corerrcor_mixedmodel,p_corerrcor_mixedmodel] = calc_ersp_linearmixedmodel(T,corerspdata_corerrcor,errcorerspdata_corerrcor,numTimes,numFreqs);

