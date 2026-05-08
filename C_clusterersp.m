%% Get average ERSP/ITC from interaction clusters 
% amb895 05/07/2026
clear;

% add paths to eeglab and limo tools
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/Resources/slanCM/')

%% Directories and names
resultsdir = hera('Abby/AS_EEG/Results');
tfcedir = hera('Abby/AS_EEG/ERSPdata/Cluster_stats');
erspdir = hera('Abby/AS_EEG/ERSPdata');
eventtype = 'prep';
basenorm = 'avgtrl';
% basenorm = 'singletrl';
frow = [4 5 6 7 37 38 39 40 41];

% correct AS ersp directory
corASerspdir = [erspdir '/corAS_data/' eventtype];
VGSerspdir = [erspdir '/VGS_data/' eventtype];

%% Load correct AS and VGS data
corASersp = load(fullfile(corASerspdir,['corAS_' eventtype '_' basenorm '_ersp.mat']));
VGSersp = load(fullfile(VGSerspdir,['VGS_' eventtype '_' basenorm '_ersp.mat']));
times = corASersp.times;
freqs  =corASersp.freqs;

corASitc = load(fullfile(corASerspdir,['corAS_' eventtype '_' basenorm '_itc.mat']));
VGSitc = load(fullfile(VGSerspdir,['VGS_' eventtype '_' basenorm '_itc.mat']));

corASid = load(fullfile(corASerspdir,['corAS_' eventtype '_' basenorm '_ids.mat']));
VGSid = load(fullfile(VGSerspdir,['VGS_' eventtype '_' basenorm '_ids.mat']));

% remove non viable from AS 
corASersp.allersp_viable = corASersp.allersp(corASid.idmat(:,4)==1 ,:,:,:);
corASitc.allitc_viable = corASitc.allitc(corASid.idmat(:,4)==1,:,:,:);
corASid.idmat_viable = corASid.idmat(corASid.idmat(:,4)==1,:);

% average ersp and itc across frow
corASersp.allersp_frow = squeeze(mean(corASersp.allersp_viable(:,:,:,frow),4));
VGSersp.allersp_frow = squeeze(mean(VGSersp.allersp(:,:,:,frow),4));
corASitc.allitc_frow = squeeze(mean(corASitc.allitc_viable(:,:,:,frow),4));
VGSitc.allitc_frow = squeeze(mean(VGSitc.allitc(:,:,:,frow),4));

% set up data to loop over
data.ersp.corAS = corASersp.allersp_frow;
data.ersp.VGS = VGSersp.allersp_frow;
data.itc.corAS = corASitc.allitc_frow;
data.itc.VGS = VGSitc.allitc_frow;

% load merged7t AS score table
load(fullfile(resultsdir,'merged7t_ASscore_table.mat'))   

% define temp table
temptable = merged7t_ASscore_table;

% total subjects
N = height(merged7t_ASscore_table);

% set up table to define time and frequency duration of each cluster
clustertable = table('Size',[0,7],'VariableTypes',["string" "string" repmat("double",1,5)],...
    'VariableNames',["modelterm" "datatype" "minfreq","maxfreq","mintime","maxtime","clusternum"]);
%% Loop over ERSP/ ITC condition, age, and interaction clusters
modelterms = ["invage", "condition", "interaction"];
datatypes = ["ersp", "itc"];

for d = 1:length(datatypes)
    currdatatype = char(datatypes(d));
    for m = 1:length(modelterms)
        currterm = modelterms(m);
        % print out progress
        fprintf('%s %s\n',currdatatype,currterm)
        % load in tfce data
        tfce = load(fullfile(tfcedir,['corASvsVGS_tfceresults_' currdatatype '_' eventtype '_' basenorm '.mat']));
        currclustermask = tfce.tfce_results.(currterm).mask;
        % get number of clusters
        [labeledclusters ,numclusters] = bwlabel(currclustermask);
        % view clusters
        figure; surf(labeledclusters,'EdgeColor','none');view(2); colorbar    
        % print out number of clusters
        fprintf('%s %s total has %d clusters\n',currdatatype,currterm,numclusters)
        
        % determine size of each cluster
        clusterstats = regionprops(labeledclusters);
        % if cluster size has an area less than 70 pixels remove cluster from analyses
        toosmall_clustnums = find(([clusterstats.Area] < 70) == 1);
        % loop over clusters to be removed
        for p = 1:length(toosmall_clustnums)
            currnum = toosmall_clustnums(p);
            % in labeledclusters set the cluster to be removed to zero
            [row,col] = find(labeledclusters == currnum);
            for i = 1:length(row)
                labeledclusters(row(i),col(i)) = 0;
            end
        end
        % update number of clusters
        newclustermask = labeledclusters > 0;
        [labeledclusters,numclusters] = bwlabel(newclustermask);

        % print out clusters removed
        fprintf("%d clusters removed from total %s\n",length(toosmall_clustnums),currterm)         
        figure; surf(labeledclusters,'EdgeColor','none');view(2); colorbar    
        
        %% loop over clusters
        for c = 1:numclusters
            % print out current cluster
            fprintf("%d\n",c)
            % current cluster mask
            currclustermask = c == labeledclusters;
            % determine frequency band and start and stop times of cluster
           [freqIDX,timeIDX] = find(currclustermask);
           C = height(clustertable) + 1;
           clustertable(C,"modelterm") = table(string(currterm));
           clustertable(C,"datatype") = table(string(currdatatype));
           clustertable(C,"minfreq") = table(freqs(min(freqIDX)));
           clustertable(C,"maxfreq") = table(freqs(max(freqIDX)));
           clustertable(C,"mintime") = table(times(min(timeIDX)));
           clustertable(C,"maxtime") = table(times(max(timeIDX)));
           clustertable(C,"clusternum") = table(c);

           % define ersp/ power vectors
           avgcorAS_data = zeros(N,1);
           avgVGS_data = zeros(N,1);

           % loop through all subjects in merged7t_ASscore_table
           for s = 1:N
               % get lunaid and eeg_date from merged7t_ASscore_table
               lunaid = merged7t_ASscore_table.lunaid(s);
               eeg_date = merged7t_ASscore_table.eeg_date(s);

               %% CorAS
               % find subject's ersp/power file (correct AS and VGS)
               corASdata_idx=find(corASid.idmat_viable(:,1) == lunaid & corASid.idmat_viable(:,2)==eeg_date);
               %% Average ERSP 
               if isempty(corASdata_idx)
                   % if no file found set power and ersp to NanN
                   avgcorAS_data(s,1) = NaN; 
               else
                   %% FIX ME instead use data already loaded above
                   sub_corASdata = squeeze(data.(currdatatype).corAS(corASdata_idx,:,:));
%                    corASerspfile = corASerspfiles(corASerspfileidx);
                   % load file
%                    substructcorASersp = load(fullfile(corASerspfile.folder,corASerspfile.name),...
%                        'times','freqs','sub_totersp');
%                    ERSPTIMES = substructcorASersp.times;

%                     % average across f-row
%                    sub_toterspcorAS = squeeze(mean(...
%                        substructcorASersp.sub_totersp(:,...
%                        :,:,[4 5 6 7 37 38 39 40 41]),4));

                   % mask ersp and power with current cluster mask and average value and save to vectors
                   avgcorAS_data(s,1) = mean(nonzeros(sub_corASdata.*currclustermask));


               end


                %% VGS
                % find subject's ersp/power file (correct AS and VGS)
               VGSdata_idx=find(VGSid.idmat(:,1) == lunaid & VGSid.idmat(:,2)==eeg_date);
               %% Average ERSP 
               if isempty(VGSdata_idx)
                   % if no file found set power and ersp to NanN
                   avgVGS_data(s,1) = NaN; 
               else
                   sub_VGSdata = squeeze(data.(currdatatype).VGS(VGSdata_idx,:,:));
                   avgVGS_data(s,1) = mean(nonzeros(sub_VGSdata.*currclustermask));
               end
                % find subject's ersp/power file (correct AS and VGS)
%                VGSerspfilenames = {VGSerspfiles.name};
%                VGSersp_isfile = contains(VGSerspfilenames,[num2str(lunaid) '_' num2str(eeg_date)]);
%                VGSerspfileidx = find(VGSersp_isfile);
%                %% Average ERSP 
%                if isempty(VGSerspfileidx)
%                    % if no file found set power and ersp to NanN
%                    avgVGS_data(s,1) = NaN; 
%                else
%                    VGSerspfile = VGSerspfiles(VGSerspfileidx);
%                    % load file
%                    substructVGSersp = load(fullfile(VGSerspfile.folder,VGSerspfile.name),...
%                        'times','freqs','sub_totersp');
%                    ERSPTIMES = substructVGSersp.times;
% 
%                     % average across f-row
%                    sub_toterspVGS = squeeze(mean(...
%                        substructVGSersp.sub_totersp(:,...
%                        :,:,[4 5 6 7 37 38 39 40 41]),4));
% 
%                    % mask ersp and power with current cluster mask and average value and save to vectors
%                    avgVGS_data(s,1) = mean(nonzeros(sub_toterspVGS.*currclustermask));
% 
%                end


           end
           % add vectors of ersp and power for corAS and VGS to big table
           powertable_temp = table(avgcorAS_data, avgVGS_data,...
               'VariableNames',...
               {sprintf('corAS_c%d_%s_%s',c,currdatatype,currterm),...
               sprintf('VGS_c%d_%s_%s',c,currdatatype,currterm)});
           % join with merged7t_ASscore_table
           temptable = [temptable,powertable_temp];    

        end    

    end
    
end

%% save table
currdate = datetime();
currdate.Format = 'yyyyMMdd';
clusterpowertable = temptable;
save(fullfile(resultsdir,['ersp_itc_clusters_',basenorm,'_',char(currdate),'.mat']),"clusterpowertable")
writetable(clusterpowertable,fullfile(resultsdir,['ersp_itc_clusters_',basenorm,'_',char(currdate),'.csv']))
save(fullfile(resultsdir,['ersp_itc_clusterinfo_' basenorm, '_' char(currdate) '.mat']),"clustertable")
writetable(clustertable,fullfile(resultsdir,['ersp_itc_clusterinfo_' basenorm '_' char(currdate) '.csv']))

