clear
% load in necessary paths and AS score data 
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/'));
addpath(genpath('/resources/Euge/'))
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable.mat')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/')
%%
maindir = hera('Abby/preprocessed_data');
task='anti';
taskdirectory = [maindir '/' task]; 
setfilesDir = [taskdirectory, '/AfterWhole/ICAWholeClean/*_icapru.set'];
filesDir = dir(setfilesDir);
% getting all files
for i = 1:length(filesDir)
    folder = filesDir(i).folder;
    name = filesDir(i).name;
    setfiles{i} = [folder,'/',name];
end
%%
n = length(setfiles);
sz = size(ErrorLatencyTable);
ICAvarAgeTable = table('Size',[sz(1),6],'VariableTypes',...
    {'double','double','double','string','double','double'},...
    'VariableNames',{'LunaID','ScanDate','Age','Sex','BeforeICAVar','AfterICAVar'});
for i=1:n
    inputfile = setfiles{i};
    [d, currentName, ext ] = fileparts(inputfile);
    parts = split(currentName,'_');
    subid = str2double(parts{1});
    scandate = str2double(parts{2});
    idx = find(ErrorLatencyTable.("Luna ID")==subid & ErrorLatencyTable.("Scan Date")==scandate);
    if isempty(idx)
        continue
    end
    EEG = pop_loadset(inputfile);
    varBefore = EEG.etc.varBeforeICArej;
    varAfter = EEG.etc.varAfterICArej;
    ageAtScan = ErrorLatencyTable.("Age")(idx);
    sex = string(table2cell(ErrorLatencyTable(idx,"Sex")));
    ICAvarAgeTable(i,:) = {subid,scandate,ageAtScan,sex,varBefore,varAfter};
end
idxMissingData = find(ICAvarAgeTable.("LunaID") == 0);
ICAvarAgeTable(idxMissingData,:) = [];

%%
merged7t = readtable('merged_7t.csv'); % to get ages
asScore = load('eeg_anti.mat');
merge7tIDs = unique(merged7t.lunaid);
errorLatIDs = unique(ErrorLatencyTable.("Luna ID"));
asScoreIDs = unique(asScore.datatable.("LunaID"));
