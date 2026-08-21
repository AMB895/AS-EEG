%% Compute ERSP stats
% amb895 05/07/2026
clear;
% Full model: ersp/ itc ~ condition * invage + (1|lunaid)
%% Directories to data, paths to add and variables to define
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
maindir = hera('Abby/AS_EEG/ERSPdata');
eventtype = 'prep';
basenorm = 'avgtrl';
% basenorm = 'singletrl';
corASdir = [maindir '/corAS_data/' eventtype '/'];
VGSdir = [maindir '/VGS_data/' eventtype '/'];
errcorASdir = [maindir '/errcorAS_data/' eventtype '/'];
% datatypes = ["ersp","itc"];
datatypes = "ersp";
row = 'prow';
% electrodes =  'frow';

RERUN = 0;
if strcmp(row,'frow')
    electrodes = [4 5 6 7 37 38 39 40 41];
elseif strcmp(row,'prow')
    electrodes = [19 20 21 22 23 31 55 56 57 58 59];
end

%% Load correct AS and VGS data
corASersp = load(fullfile(corASdir,['corAS_' eventtype '_' basenorm '_ersp.mat']));
VGSersp = load(fullfile(VGSdir,['VGS_' eventtype '_' basenorm '_ersp.mat']));
times = corASersp.times;
freqs  =corASersp.freqs;

corASitc = load(fullfile(corASdir,['corAS_' eventtype '_' basenorm '_itc.mat']));
VGSitc = load(fullfile(VGSdir,['VGS_' eventtype '_' basenorm '_itc.mat']));

corASid = load(fullfile(corASdir,['corAS_' eventtype '_' basenorm '_ids.mat']));
VGSid = load(fullfile(VGSdir,['VGS_' eventtype '_' basenorm '_ids.mat']));

% remove non viable from AS 
corASersp.allersp_viable = corASersp.allersp(corASid.idmat(:,4)==1 ,:,:,:);
corASitc.allitc_viable = corASitc.allitc(corASid.idmat(:,4)==1,:,:,:);
corASid.idmat_viable = corASid.idmat(corASid.idmat(:,4)==1,:);

% average ersp and itc across row
corASersp.allersp_row = squeeze(mean(corASersp.allersp_viable(:,:,:,electrodes),4));
VGSersp.allersp_row = squeeze(mean(VGSersp.allersp(:,:,:,electrodes),4));
corASitc.allitc_row = squeeze(mean(corASitc.allitc_viable(:,:,:,electrodes),4));
VGSitc.allitc_row = squeeze(mean(VGSitc.allitc(:,:,:,electrodes),4));

% variables to loop over
data.ersp.corAS = corASersp.allersp_row;
data.ersp.VGS = VGSersp.allersp_row;
data.itc.corAS = corASitc.allitc_row;
data.itc.VGS = VGSitc.allitc_row;

%% loop over data types (ersp and itc)
for p = 1:length(datatypes)
    currdatatype = char(datatypes(p));
    % skip if created
    if exist(fullfile(maindir,'/Cluster_stats/',['corASvsVGS_' char(currdatatype) '_' eventtype '_' basenorm '_' row '.mat']),'file') && RERUN==0
        fprintf('Computed stats for %s %s %s; skipping\n',currdatatype, eventtype, basenorm)
    else
        % print out ersp or itc
        fprintf("Computing stats for %s %s %s\n",currdatatype,eventtype,basenorm)
        % get data type
        corASdata = data.(currdatatype).corAS;
        VGSdata = data.(currdatatype).VGS;
        % set up tables
        varnames = ["lunaid" "invage" "condition" "data"];
        vartypes = ["categorical" "double" "categorical" "double"];
        % correct AS table
        Tcor = table('Size',[length(corASid.idmat_viable) 4],'VariableTypes',vartypes,'VariableNames',varnames);
        Tcor.lunaid = corASid.idmat_viable(:,1);
        Tcor.invage = 1./corASid.idmat_viable(:,3);
        Tcor(:,"condition") = table(categorical(cellstr('corAS')));
        % VGS table
        Tvgs = table('Size',[length(VGSid.idmat) 4],'VariableTypes',vartypes,'VariableNames',varnames);
        Tvgs.lunaid = VGSid.idmat(:,1);
        Tvgs.invage = 1./VGSid.idmat(:,3);
        Tvgs(:,"condition") = table(categorical(cellstr('VGS')));
 
        % loop over time points
        for t = 1:length(times)
            fprintf("%.1f ms\n",times(t))
            % loop over frequency points
            for f = 1:length(freqs)
                Tcor.data = squeeze(corASdata(:,f,t));
                Tvgs.data = squeeze(VGSdata(:,f,t));
                % combine tables
                T = [Tvgs;Tcor];
                % run model at time/freq point
                templme = fitlme(T, 'data ~ invage*condition + (1|lunaid)','FitMethod','REML');
                coeffs = table2cell(dataset2table(templme.Coefficients));
                % Get t-valeus
                tval.intercept(f,t) = coeffs{1,4};
                tval.invage(f,t) = coeffs{2,4};
                tval.condition(f,t) = coeffs{3,4};
                tval.interaction(f,t) = coeffs{4,4};                
            end
        end
        % save t-values
        save(fullfile(maindir,'/Cluster_stats/',['corASvsVGS_' char(currdatatype) '_' eventtype '_' basenorm '_' row '.mat']),...
            'tval','times','freqs')
    end
end

