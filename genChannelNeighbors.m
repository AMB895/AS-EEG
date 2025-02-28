%% Threshold-Free Cluster Enhancement Method
% From limo_tools
% Functions:
% [tfce_score, thresholded_maps] = limo_tfce(type, data,channeighbstructmat, updatebar, E, H, dh)
%   type: 1 for 1D data, 2 for 2D data, 3 for 3D data
%   data: map of t/F values or set of t/F maps computed under H0
%   channeighbstructmat: neighbourhood matrix for clustering (size is # of chans x # of chans)

% adding necessary paths
addpath('/Volumes/Hera/Abby/AS_EEG/STUDY/') % to get ERSP,ITC & ERP data
addpath('/Volumes/Hera/Abby/AS_EEG/limo_tools/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/preprocessed_data/anti/AfterWhole/epochclean_homogenize/')

eeglab;

% loading data
load('erspSTUDY.mat')
load('itcSTUDY.mat')
erspdata = erspdata{1,1};
itcdata = itcdata{1,1};

% STUDY file path and name
studypath = '/Volumes/Hera/Abby/AS_EEG/STUDY/';
studyname = 'AS_correct_PrepPeriod.study';

% % load study
 [STUDY, ALLEEG] = pop_loadstudy('filename',studyname,'filepath',studypath);

%% get channel neighbours
[~,hostname] = system('hostname');
beginningHostname = hostname(1:4);
if ~strcmp(beginningHostname,'rhea') % MEX error need to run on rhea
    fprintf('Need to run on Rhea\n')
else
    % getting channel neighbours
    [STUDY] = pop_statparams(STUDY,'groupstats','on','condstats','on','mode','fieldtrip','fieldtripmcorrect','cluster');
    [STUDY ,neighbours] = std_prepare_neighbors(STUDY,ALLEEG,'force','on');
    [STUDY, ALLEEG] = pop_savestudy(STUDY,[],'filename',studyname,'filepath',studypath);
    save(fullfile(studypath,'neighbourStruct.mat'),'neighbours')
end

%% created channel neighbour logical table and matrix
% load structure and generate channel neighbour logical matrix
load(fullfile(studypath,'neighbourStruct.mat'))
numChans = length(neighbours);
channeighbmatrix = zeros(numChans);
varTypes = repmat({'double'},1,numChans);
varNames = {neighbours.label};
channeighbtable = table('Size',size(channeighbmatrix),'VariableTypes',varTypes,'VariableNames',varNames,...
    'RowNames',varNames);
for chan = 1:numChans
    chanName = neighbours(chan).label;
    chanNeighbs = neighbours(chan).neighblabel;
    numNeighbs = length(chanNeighbs);
    for neighb = 1:numNeighbs
        neighbName = chanNeighbs(neighb);
        channeighbtable(neighbName,chan) = {1};
    end
end
channeighbmatrix = cell2mat(table2cell(channeighbtable));
save(fullfile(studypath,'logicalchanneighbmatrix.mat'),'channeighbmatrix')

