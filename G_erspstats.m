%% Compute ERSP stats
% amb895 05/07/2026
clear;
% Full model: ersp/ itc ~ condition * invage + (1|lunaid)
%% Directories to data, paths to add and variables to define
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
maindir = hera('Abby/AS_EEG/ERSPdata');
eventtype = 'prep';
% basenorm = 'avgtrl';
basenorm = 'singletrl';
corASdir = [maindir '/corAS_data/' eventtype '/'];
VGSdir = [maindir '/VGS_data/' eventtype '/'];
errcorASdir = [maindir '/errcorAS_data/' eventtype '/'];
datatypes = ["ersp","itc"];
RERUN = 0;
frow = [4 5 6 7 37 38 39 40 41];
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

% average ersp and itc across frow
corASersp.allersp_frow = squeeze(mean(corASersp.allersp_viable(:,:,:,frow),4));
VGSersp.allersp_frow = squeeze(mean(VGSersp.allersp(:,:,:,frow),4));
corASitc.allitc_frow = squeeze(mean(corASitc.allitc_viable(:,:,:,frow),4));
VGSitc.allitc_frow = squeeze(mean(VGSitc.allitc(:,:,:,frow),4));

% variables to loop over
data.ersp.corAS = corASersp.allersp_frow;
data.ersp.VGS = VGSersp.allersp_frow;
data.itc.corAS = corASitc.allitc_frow;
data.itc.VGS = VGSitc.allitc_frow;

%% loop over data types (ersp and itc)
for p = 1:length(datatypes)
    currdatatype = char(datatypes(p));
    % skip if created
    if exist(fullfile(maindir,'/Cluster_stats/',['corASvsVGS_' char(currdatatype) '_' eventtype '_' basenorm '.mat']),'file') && RERUN==0
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
        save(fullfile(maindir,'/Cluster_stats/',['corASvsVGS_' char(currdatatype) '_' eventtype '_' basenorm '.mat']),...
            'tval','times','freqs')
    end
end

