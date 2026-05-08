%% Run TFCE on ersp/itc stats
% amb895 05/07/2026

% add paths to eeglab and limo tools
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/Resources/slanCM/') 

% define variables
datatypes = ["ersp" "itc"];
eventtype = 'prep';
basenorm = 'avgtrl';

% define directories and file names
maindir = hera('Abby/AS_EEG/ERSPdata/Cluster_stats');

% rerun= 0; no rerun / rerun = 1; rerun
RERUN =0;

% plot=0; no plot / plot=1; plot
PLOT=1;

%% loop over data types
% Model 1: ERSP/ITC ~ condition*inverse age 
for p = 1:length(datatypes)
   currdatatype = char(datatypes(p));
   if exist(fullfile(maindir,['corASvsVGS_tfceresults_' currdatatype '_' eventtype '_' basenorm '.mat']),'file') && RERUN==0
        fprintf('Computed TFCE for Full Model %s %s %s; skipping\n',currdatatype,eventtype,basenorm)
        if PLOT
            load(fullfile(maindir,['corASvsVGS_tfceresults_' currdatatype '_' eventtype '_' basenorm '.mat']))
            load(fullfile(maindir, ['corASvsVGS_' currdatatype '_' eventtype '_' basenorm '.mat']))
        end
   else
       % print progress
       fprintf("%s %s\n", currdatatype,basenorm)
       % load current data stats
       filename = ['corASvsVGS_' currdatatype '_' eventtype '_' basenorm '.mat'];
       % structure: tval, times,freqs
       load(fullfile(maindir,filename))
       
       % TFCE on t-values
        tfce_results.invage.tfcescores = limo_tfce(2, tval.invage,[]);
        tfce_results.condition.tfcescores = limo_tfce(2, tval.condition, []); 
        tfce_results.interaction.tfcescores = limo_tfce(2, tval.interaction, []);
        
        % Permute t-values and run TFCE
        tfce_results.invage.permscores = calc_perm_tfce_2d(tval.invage,1000); 
        tfce_results.condition.permscores = calc_perm_tfce_2d(tval.condition,1000);
        tfce_results.interaction.permscores = calc_perm_tfce_2d(tval.interaction,1000);

        % Calculate p<0.01 mask
        tfce_results.invage.mask = calc_tfce_mask_2d(tfce_results.invage.tfcescores, tfce_results.invage.permscores,0.01);
        tfce_results.condition.mask = calc_tfce_mask_2d(tfce_results.condition.tfcescores, tfce_results.condition.permscores, 0.01);
        tfce_results.interaction.mask = calc_tfce_mask_2d(tfce_results.interaction.tfcescores, tfce_results.interaction.permscores, 0.01);      

        % save tfce results
        save(fullfile(maindir,['corASvsVGS_tfceresults_' currdatatype '_' eventtype '_' basenorm '.mat']), ...
            "tfce_results", "times","freqs")
        
   end
   
   %% Plotting
    if PLOT
        
        % define MAP
        MAP = slanCM('viridis',256);

        % plot clusters with t-values
        % interaction clusters
        mask_interaction = tfce_results.interaction.mask;
            plot_tfce_clusters(times,freqs,tval.interaction,'mask',mask_interaction,'title',['Interaction Effect ' currdatatype],'subtitle','F-row electrodes',...
                'xlab','Preparatory Time (ms)','ylab','Frequency (Hz)','clab','t-value','addlines',1,'alpha',0.6,'map',MAP)
            
        mask_condition = tfce_results.condition.mask;
        plot_tfce_clusters(times,freqs,tval.condition,'mask',mask_condition,'title',['Condition Effect ' currdatatype],'subtitle','F-row electrodes',...
                'xlab','Preparatory Time (ms)','ylab','Frequency (Hz)','clab','t-value','addlines',1,'alpha',0.6,'map',MAP)
        
        mask_age = tfce_results.invage.mask;
        plot_tfce_clusters(times,freqs,-1.*tval.invage,'mask',mask_age,'title',['Age Effect ' currdatatype],'subtitle','F-row electrodes',...
                'xlab','Preparatory Time (ms)','ylab','Frequency (Hz)','clab','t-value','addlines',1,'alpha',0.6,'map',MAP)
         
    end

   
end

%% Functions
function [max_perm_tfce] = calc_perm_tfce_2d(tvals,nperm)
size1 = size(tvals,1);
size2 = size(tvals,2);
max_perm_tfce = zeros(nperm,1);
parfor n=1:nperm
    % getting total size of t matrix
    totalSize = size1*size2;
    
    % reshape t matrix
    reshape_tvals = reshape(tvals,[1,totalSize]);
    
    % permute t matrix
    permtvals = reshape_tvals(randperm(totalSize)); 
    
    % reshape back to original size
    permtvals = reshape(permtvals,[size1,size2]);
        
    % run limo_tfce on permuted t matrix
    perm_tfce_scores = limo_tfce(2,permtvals,[],0);
    
    % save max permutation tfce score per iteration
    max_perm_tfce(n,:) = max(perm_tfce_scores,[],[1 2]);
end

end

function [mask] = calc_tfce_mask_2d(tfcescores,permtfcescores,p_thres)
% determine pth percentile
percentile = 1-p_thres;

% determine pth percentile for permuted TFCE scores
n = length(permtfcescores)*percentile;

% sort permuted TFCE scores in acending order
sorted_permTFCEscores = sort(permtfcescores);

% Threshold for 95th percentile
thres = sorted_permTFCEscores(n);

% Threshold Mask
mask = tfcescores > thres;
end