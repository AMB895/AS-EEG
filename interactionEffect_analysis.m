%% Average power w/in interaction effect (Age:Trial Type) for Correct vs. VGS
% add paths
addpath('/Volumes/Hera/Abby/AS_EEG/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ClusterStats/')

%% Loading in ERSP data and significance mask for interaction effect (Correct  AS vs. VGS)
clear; close all;
% load error latency table
load('ErrorLatencyTable_20250320.mat')

% load correct ersp data and ID mat
load('corERSPdata.mat')
load('corIDmatrix.mat')

% remove participants that are not viable (Correct AS)
numSubs = size(corIDmatrix,1);

corerspdata_viable=[];
corpowbase_viable=[];
corIDmatrix_viable=[];
for currentSub = 1:numSubs
    id = corIDmatrix(currentSub,1);
    scandate = corIDmatrix(currentSub,2);
    idx = find(ErrorLatencyTable.LunaID == id & ErrorLatencyTable.ScanDate==scandate);
    subtable = ErrorLatencyTable(idx,:);
    isviable = subtable.Viable==1;
    if isviable
        corerspdata_viable(end+1,:,:,:) = corerspdata(currentSub,:,:,:);
        corpowbase_viable(end+1,:,:) = corpowbase(currentSub,:,:);
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
corerspdata_viable_frow = squeeze(mean(corerspdata_viable(:,[4 5 6 7 37 38 39 40 41],:,:),2));
corpowbase_viable_frow = squeeze(mean(corpowbase_viable(:,[4 5 6 7 37 38 39 40 41],:),2));

% load in vgs ersp data and ID matrix
load('vgsERSPdata.mat')
load('vgsIDmatrix.mat')

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
vgserspdata_frow = squeeze(mean(vgserspdata(:,[4 5 6 7 37 38 39 40 41],:,:),2));
vgspowbase_frow = squeeze(mean(vgspowbase(:,[4 5 6 7 37 38 39 40 41],:),2));

% Match IDs from correct AS and VGS trials
corvgsIDmatrix = [];
vgserspdata_corvgs = [];
corerspdata_corvgs = [];
vgspowbase_corvgs = [];
corpowbase_corvgs = [];
for idxVGS = 1:size(vgsIDmatrix,1)
   currentSub = vgsIDmatrix(idxVGS,1);
   currentDate = vgsIDmatrix(idxVGS,2);
   idxAS = find(corIDmatrix_viable(:,1)==currentSub & corIDmatrix_viable(:,2) == currentDate);
   if ~isempty(idxAS) % there is a matching ID and date
       corvgsIDmatrix(end+1,:) = vgsIDmatrix(idxVGS,:);
       vgserspdata_corvgs(end+1,:,:) = vgserspdata_frow(idxVGS,:,:);
       corerspdata_corvgs(end+1,:,:) = corerspdata_viable_frow(idxAS,:,:);
       vgspowbase_corvgs(end+1,:) = vgspowbase_frow(idxVGS,:);
       corpowbase_corvgs(end+1,:) = corpowbase_viable_frow(idxAS,:);
   end
end

% load in mask
load('corvgsMixedModelClusters.mat')

%% generating masks for each cluster
% number the clusters
interaction_clusters_map_corvgs = mask_corvgs_mixedmodel.interaction.*t_corvgs_mixedmodel.interaction;
[labeled_clusters, numClusters] = bwlabel(interaction_clusters_map_corvgs);
numFreqs = length(freqs);
numTimes = length(times);
numSubs = size(corvgsIDmatrix,1);

% setting up matricies to convert to table
Tdata_pow = zeros(numSubs,11);
Tdata_pow(:,1) = corvgsIDmatrix(:,1);
Tdata_pow(:,2) = corvgsIDmatrix(:,4);
Tdata_pow(:,3) = corvgsIDmatrix(:,3);

Tdata_sigchange = zeros(numSubs,11);
Tdata_sigchange(:,1) = corvgsIDmatrix(:,1);
Tdata_sigchange(:,2) = corvgsIDmatrix(:,4);
Tdata_sigchange(:,3) = corvgsIDmatrix(:,3);

% Table set up:
% [ID, VisitNum, Age, C1 Pow, C2 Pow,C3 Pow, C4 Pow,
% C1 Signal Change, C2 Signal Change,C3 Signal Change,C4 Signal change]

% generate mask for each cluster
for currentCluster = 1:numClusters
    % find individual clusters
    cluster_mask_corvgs = labeled_clusters == currentCluster;
    % determine frequency band
    [freqIdx,timeIdx] = find(cluster_mask_corvgs);
    minFreqIdx = min(freqIdx);
    maxFreqIdx = max(freqIdx);
    minFreq = freqs(minFreqIdx);
    maxFreq = freqs(maxFreqIdx);
    
     % plot individual clusters
    figure;
    surf(times,freqs,double(cluster_mask_corvgs),'EdgeColor','none')
    view(2);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)')
    title(sprintf('Cluster %d %2.1f-%2.1f Hz',currentCluster,minFreq,maxFreq))
    
    % save frequency bands
    clusterFreqbands(currentCluster,1) = minFreq;
    clusterFreqbands(currentCluster,2) = maxFreq;
    
    % average the baseline power within the frequency band
    avgcorpowbase_corvgs = mean(corpowbase_corvgs(:,minFreqIdx:maxFreqIdx),2);
    avgvgspowbase_corvgs = mean(vgspowbase_corvgs(:,minFreqIdx:maxFreqIdx),2);
    
    for currentSub = 1:numSubs
        % get current subject's baseline power and ersp data
        subavgcorpowbase = avgcorpowbase_corvgs(currentSub);
        subavgvgspowbase = avgvgspowbase_corvgs(currentSub);
        subcorersp = squeeze(corerspdata_corvgs(currentSub,:,:));
        subvgsersp = squeeze(vgserspdata_corvgs(currentSub,:,:));
        
        % mask the ersp data with current cluster mask
        subcorersp_cluster = cluster_mask_corvgs.*subcorersp;
        subvgsersp_cluster = cluster_mask_corvgs.*subvgsersp;
        
        % get overall power of ERSP image
        overallcorpow = mean(subcorersp,'all');
        overallvgspow = mean(subvgsersp,'all');
        
        % average power within the cluster (total cluster power)
        totcorclusterpow = mean(nonzeros(subcorersp_cluster));
        totvgsclusterpow = mean(nonzeros(subvgsersp_cluster));
        
        % relative cluster power
        relcorclusterpow = totcorclusterpow/overallcorpow;
        relvgsclusterpow = totvgsclusterpow/overallvgspow;
        
        % add total cluster power to matrix
        Tdata_pow(currentSub,3+currentCluster) = totcorclusterpow;
        Tdata_pow(currentSub,11+currentCluster) = totvgsclusterpow;
        
        % add relative cluster power to matrix
        Tdata_pow(currentSub,7+currentCluster) = relcorclusterpow;
        Tdata_pow(currentSub,15+currentCluster) = relvgsclusterpow;
        
        % calculate % signal change from baseline
        corsigchange = totcorclusterpow/subavgcorpowbase;
        vgssigchange = totvgsclusterpow/subavgvgspowbase;
        
        % add % signal change to matrix
        Tdata_sigchange(currentSub,3+currentCluster) = corsigchange;
        Tdata_sigchange(currentSub,3+numClusters+currentCluster) = vgssigchange;
    end
end

% Table set up:
% [ID, VisitNum, Age, C1 Pow, C2 Pow,C3 Pow, C4 Pow,
% C1 Signal Change, C2 Signal Change,C3 Signal Change,C4 Signal change]

% setting up tables
varNames = {'ID','VisitNum','Age','totcorC1Pow','totcorC2Pow','totcorC3Pow','totcorC4Pow',...
    'relcorC1Pow','relcorC2Pow','relcorC3Pow','relcorC4Pow',...
    'totvgsC1Pow','totvgsC2Pow','totvgsC3Pow','totvgsC4Pow',...
    'relvgsC1Pow','relvgsC2Pow','relvgsC3Pow','relvgsC4Pow'};
T_pow = table('Size',[numSubs 3+4*numClusters],'VariableTypes',repmat("double",[1 3+4*numClusters]),...
    'VariableNames',varNames);
T_pow.ID = Tdata_pow(:,1); 
T_pow.VisitNum = Tdata_pow(:,2);
T_pow.Age = Tdata_pow(:,3);
% add correct total and relative power to table
T_pow.totcorC1Pow = Tdata_pow(:,4); 
T_pow.totcorC2Pow = Tdata_pow(:,5); 
T_pow.totcorC3Pow = Tdata_pow(:,6);
T_pow.totcorC4Pow = Tdata_pow(:,7);
T_pow.relcorC1Pow = Tdata_pow(:,8); 
T_pow.relcorC2Pow = Tdata_pow(:,9); 
T_pow.relcorC3Pow = Tdata_pow(:,10);
T_pow.relcorC4Pow = Tdata_pow(:,11);

% add vgs total and relative power to table
T_pow.totvgsC1Pow = Tdata_pow(:,12);
T_pow.totvgsC2Pow = Tdata_pow(:,13);
T_pow.totvgsC3Pow = Tdata_pow(:,14);
T_pow.totvgsC4Pow = Tdata_pow(:,15);

T_pow.relvgsC1Pow = Tdata_pow(:,16);
T_pow.relvgsC2Pow = Tdata_pow(:,17);
T_pow.relvgsC3Pow = Tdata_pow(:,18);
T_pow.relvgsC4Pow = Tdata_pow(:,19);

varNames = {'ID','VisitNum','Age','corC1SigChange','corC2SigChange','corC3SigChange','corC4SigChange','vgsC1SigChange','vgsC2SigChange','vgsC3SigChange','vgsC4SigChange'};
T_sigchange = table('Size',[numSubs 3+2*numClusters],'VariableTypes',repmat("double",[1 3+2*numClusters]),...
    'VariableNames',varNames);
T_sigchange.ID = Tdata_sigchange(:,1);
T_sigchange.VisitNum = Tdata_sigchange(:,2); 
T_sigchange.Age = Tdata_sigchange(:,3);
T_sigchange.corC1SigChange = Tdata_sigchange(:,4); 
T_sigchange.corC2SigChange = Tdata_sigchange(:,5); 
T_sigchange.corC3SigChange = Tdata_sigchange(:,6);
T_sigchange.corC4SigChange = Tdata_sigchange(:,7);
T_sigchange.vgsC1SigChange = Tdata_sigchange(:,8); 
T_sigchange.vgsC2SigChange = Tdata_sigchange(:,9);
T_sigchange.vgsC3SigChange = Tdata_sigchange(:,10);
T_sigchange.vgsC4SigChange = Tdata_sigchange(:,11);

% save tables
writetable(T_pow,'/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/clusterpower_corvgsinteraction.csv')
writetable(T_sigchange,'/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/signalchange_corvgsinteraction.csv')


%% Average power w/in interaction effect (Age:Trial Type) for Correct vs. VGS
% add paths
addpath('/Volumes/Hera/Abby/AS_EEG/')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/vgs')
addpath('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/')

%% Loading in ERSP data and significance mask for interaction effect (Correct  AS vs. Error AS)

% load correct ersp data and ID mat
load('errcorERSPdata.mat')
load('errcorIDmatrix.mat')

% remove participants that are not viable (Error Correct AS)
numSubs = size(errcorIDmatrix,1);

errcorerspdata_viable=[];
errcorpowbase_viable=[];
errcorIDmatrix_viable=[];
for currentSub = 1:numSubs
    id = errcorIDmatrix(currentSub,1);
    scandate = errcorIDmatrix(currentSub,2);
    idx = find(ErrorLatencyTable.LunaID == id & ErrorLatencyTable.ScanDate==scandate);
    subtable = ErrorLatencyTable(idx,:);
    isviable = subtable.Viable==1;
    if isviable
        errcorerspdata_viable(end+1,:,:,:) = errcorerspdata(currentSub,:,:,:);
        errcorpowbase_viable(end+1,:,:) = errcorpowbase(currentSub,:,:);
        errcorIDmatrix_viable(end+1,:) = errcorIDmatrix(currentSub,:);
    end
end

% Add visit number to corIDmatrix
errcorIDs = unique(errcorIDmatrix_viable(:,1));
for currentSub = 1:length(errcorIDs)
    subIdx = find(errcorIDmatrix_viable(:,1)==errcorIDs(currentSub));
    subinfo = errcorIDmatrix_viable(subIdx,:);
    numVisits = size(subinfo,1);
    for currentVisit=1:numVisits
       errcorIDmatrix_viable(subIdx(currentVisit),4) = currentVisit; 
    end
end

% Average across F-row
errcorerspdata_viable_frow = squeeze(mean(errcorerspdata_viable(:,[4 5 6 7 37 38 39 40 41],:,:),2));
errcorpowbase_viable_frow = squeeze(mean(errcorpowbase_viable(:,[4 5 6 7 37 38 39 40 41],:),2));

% Match IDs from correct AS and error AS trials
corerrcorIDmatrix = [];
errcorerspdata_corerrcor = [];
corerspdata_corerrcor = [];
errcorpowbase_corerrcor = [];
corpowbase_corerrcor = [];
for idxERRCOR = 1:size(errcorIDmatrix_viable,1)
   currentSub = errcorIDmatrix_viable(idxERRCOR,1);
   currentDate = errcorIDmatrix_viable(idxERRCOR,2);
   idxCOR = find(corIDmatrix_viable(:,1)==currentSub & corIDmatrix_viable(:,2) == currentDate);
   if ~isempty(idxCOR) % there is a matching ID and date
       corerrcorIDmatrix(end+1,:) = errcorIDmatrix_viable(idxERRCOR,:);
       errcorerspdata_corerrcor(end+1,:,:) = errcorerspdata_viable_frow(idxERRCOR,:,:);
       corerspdata_corerrcor(end+1,:,:) = corerspdata_viable_frow(idxCOR,:,:);
       errcorpowbase_corerrcor(end+1,:) = errcorpowbase_viable_frow(idxERRCOR,:);
       corpowbase_corerrcor(end+1,:) = corpowbase_viable_frow(idxCOR,:);
   end
end

% load in mask
load('corerrcorMixedModelClusters.mat')

%% generating masks for each cluster
% number the clusters
interaction_clusters_map_corerrcor = mask_corerrcor_mixedmodel.interaction.*t_corerrcor_mixedmodel.interaction;
[labeled_clusters, numClusters] = bwlabel(interaction_clusters_map_corerrcor);
numFreqs = length(freqs);
numTimes = length(times);
numSubs = size(corerrcorIDmatrix,1);

% setting up matricies to convert to table
Tdata_pow = zeros(numSubs,13);
Tdata_pow(:,1) = corerrcorIDmatrix(:,1);
Tdata_pow(:,2) = corerrcorIDmatrix(:,4);
Tdata_pow(:,3) = corerrcorIDmatrix(:,3);

Tdata_sigchange = zeros(numSubs,13);
Tdata_sigchange(:,1) = corerrcorIDmatrix(:,1);
Tdata_sigchange(:,2) = corerrcorIDmatrix(:,4);
Tdata_sigchange(:,3) = corerrcorIDmatrix(:,3);

% generate mask for each cluster
for currentCluster = 1:numClusters
    cluster_mask_corerrcor = labeled_clusters == currentCluster;
    % determine frequency band
    [freqIdx,timeIdx] = find(cluster_mask_corerrcor);
    minFreqIdx = min(freqIdx);
    maxFreqIdx = max(freqIdx);
    minFreq = freqs(minFreqIdx);
    maxFreq = freqs(maxFreqIdx);
    
    % save frequency bands
    clusterFreqbands(currentCluster,1) = minFreq;
    clusterFreqbands(currentCluster,2) = maxFreq;
    
     % plot individual clusters
    figure;
    surf(times,freqs,double(cluster_mask_corerrcor),'EdgeColor','none')
    view(2);
    xlabel('Time (ms)'); ylabel('Frequency (Hz)')
    title(sprintf('Cluster %d %2.1f-%2.1f Hz',currentCluster,minFreq,maxFreq))
    
    % average the baseline power within the frequency band
    avgcorpowbase_corerrcor = mean(corpowbase_corerrcor(:,minFreqIdx:maxFreqIdx),2);
    avgerrcorpowbase_corerrcor = mean(errcorpowbase_corerrcor(:,minFreqIdx:maxFreqIdx),2);
    
    for currentSub = 1:numSubs
        % get current subject's baseline power and ersp data
        subavgcorpowbase = avgcorpowbase_corerrcor(currentSub);
        subavgerrcorpowbase = avgerrcorpowbase_corerrcor(currentSub);
        subcorersp = squeeze(corerspdata_corerrcor(currentSub,:,:));
        suberrcorersp = squeeze(errcorerspdata_corerrcor(currentSub,:,:));
        
        % mask the ersp data with current cluster mask
        subcorersp_cluster = cluster_mask_corerrcor.*subcorersp;
        suberrcorersp_cluster = cluster_mask_corerrcor.*suberrcorersp;
        
        % average power within the cluster
        totcorclusterpow = mean(nonzeros(subcorersp_cluster));
        avgerrcorclusterpow = mean(nonzeros(suberrcorersp_cluster));
        
        % add average power to matrix
        Tdata_pow(currentSub,3+currentCluster) = totcorclusterpow;
        Tdata_pow(currentSub,3+numClusters+currentCluster) = avgerrcorclusterpow;
        
        % calculate % signal change from baseline
        corsigchange = totcorclusterpow/subavgcorpowbase;
        errcorsigchange = avgerrcorclusterpow/subavgerrcorpowbase;
        
        % add % signal change to matrix
        Tdata_sigchange(currentSub,3+currentCluster) = corsigchange;
        Tdata_sigchange(currentSub,3+numClusters+currentCluster) = errcorsigchange;
    end
end

% Table set up:
% [ID, VisitNum, Age, C1 Pow, C2 Pow,C3 Pow, C4 Pow,
% C1 Signal Change, C2 Signal Change,C3 Signal Change,C4 Signal change]

% setting up tables
varNames = {'ID','VisitNum','Age','corC1Pow','corC2Pow','corC3Pow','corC4Pow','corC5Pow','errcorC1Pow','errcorC2Pow','errcorC3Pow','errcorC4Pow','errcorC5Pow'};
T_pow = table('Size',[numSubs 3+2*numClusters],'VariableTypes',repmat("double",[1 3+2*numClusters]),...
    'VariableNames',varNames);
T_pow.ID = Tdata_pow(:,1); 
T_pow.VisitNum = Tdata_pow(:,2);
T_pow.Age = Tdata_pow(:,3);
T_pow.corC1Pow = Tdata_pow(:,4); 
T_pow.corC2Pow = Tdata_pow(:,5); 
T_pow.corC3Pow = Tdata_pow(:,6);
T_pow.corC4Pow = Tdata_pow(:,7);
T_pow.corC5Pow = Tdata_pow(:,8);
T_pow.errcorC1Pow = Tdata_pow(:,9);
T_pow.errcorC2Pow = Tdata_pow(:,10);
T_pow.errcorC3Pow = Tdata_pow(:,11);
T_pow.errcorC4Pow = Tdata_pow(:,12);
T_pow.errcorC5Pow = Tdata_pow(:,13);

varNames = {'ID','VisitNum','Age','corC1SigChange','corC2SigChange','corC3SigChange','corC4SigChange','corC5SigChange','errcorC1SigChange','errcorC2SigChange','errcorC3SigChange','errcorC4SigChange','errcorC5SigChange'};
T_sigchange = table('Size',[numSubs 3+2*numClusters],'VariableTypes',repmat("double",[1 3+2*numClusters]),...
    'VariableNames',varNames);
T_sigchange.ID = Tdata_sigchange(:,1);
T_sigchange.VisitNum = Tdata_sigchange(:,2); 
T_sigchange.Age = Tdata_sigchange(:,3);
T_sigchange.corC1SigChange = Tdata_sigchange(:,4); 
T_sigchange.corC2SigChange = Tdata_sigchange(:,5); 
T_sigchange.corC3SigChange = Tdata_sigchange(:,6);
T_sigchange.corC4SigChange = Tdata_sigchange(:,7);
T_sigchange.corC5SigChange = Tdata_sigchange(:,8);
T_sigchange.errcorC1SigChange = Tdata_sigchange(:,9); 
T_sigchange.errcorC2SigChange = Tdata_sigchange(:,10);
T_sigchange.errcorC3SigChange = Tdata_sigchange(:,11);
T_sigchange.errcorC4SigChange = Tdata_sigchange(:,12);
T_sigchange.errcorC5SigChange = Tdata_sigchange(:,13);

% save tables
writetable(T_pow,'/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/clusterpower_corerrcorinteraction.csv')
writetable(T_sigchange,'/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/signalchange_corerrcorinteraction.csv')
