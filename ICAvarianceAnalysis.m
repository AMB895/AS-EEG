clear
% load in necessary paths and AS score data 
addpath(genpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/'));
addpath(genpath('/resources/Euge/'))
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
load('/Volumes/Hera/Abby/AS_EEG/ErrorLatencyTable_20250218.mat')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal/trial_data/anti/')

eeglab;
%%
maindir = hera('Abby/preprocessed_data');
task='anti';
taskdirectory = [maindir '/' task]; 
setfilesDir = [taskdirectory, '/AfterWhole/ICAWholeClean_homogenize/*_icapru.set'];
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
noASscore = [];
ICAVarMat = [];

for i=1:n
    inputfile = setfiles{i};
    [~, currentName, ~ ] = fileparts(inputfile);
    parts = split(currentName,'_');
    subid = str2double(parts{1});
    scandate = str2double(parts{2});
    idx = find(ErrorLatencyTable.("LunaID")==subid & ErrorLatencyTable.("ScanDate")==scandate);
    if isempty(idx)
        noASscore(end+1,:) = [double(subid),double(scandate)];
        continue
    end
    EEG = pop_loadset(inputfile);
    varBefore = EEG.etc.varBeforeICArej;
    varAfter = EEG.etc.varAfterICArej;
    ageAtScan = EEG.age;
    ICAVarMat(end+1,:) = [subid,scandate,ageAtScan,varBefore,varAfter];
end

ICAVarTable = table(ICAVarMat(:,1),ICAVarMat(:,2),ICAVarMat(:,3),ICAVarMat(:,4),ICAVarMat(:,5),...
    'VariableNames',{'LunaID','ScanDate','Age','VarBeforeICA','VarAfterICA'});
save('/Volumes/Hera/Abby/AS_EEG/icaVariance.mat',"ICAVarTable")
writetable(ICAVarTable,'/Volumes/Hera/Abby/AS_EEG/icaVariance.csv')