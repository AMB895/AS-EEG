clear
addpath('/Volumes/Hera/Abby/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
data = readtable('7t_eegAS_ErrorLatency.csv');
ids = unique(data.LunaID);
for i=1:length(ids)
    firstVisitIndex(i) = find(data.LunaID == ids(i),1,'first');
end
firstVisitIndex = firstVisitIndex';
dataFirstVisit = data(firstVisitIndex,:);
% Age categories:
% <13 = children (3), >=13 & <18 = adolescents (2), >=18 = adults (1)
statsData = [dataFirstVisit.LunaID,dataFirstVisit.Age,dataFirstVisit.ErrorRate_noDropped_,dataFirstVisit.AverageLatency_correctTrials_]; % Luna ID, age, sex, error rate, latency
for j = 1:length(statsData)
    age = statsData(j,2);
    if age < 13
        ageCategory(j) = 3;
    elseif age >= 13 && age <18
        ageCategory(j) = 2;
    elseif age >= 18
        ageCategory(j) = 1;
    end
end
ageCategory = ageCategory';
statsTable = table('Size',[170 6],'VariableTypes',{'double','double','string','double','double','double'},...
    'VariableNames',{'LunaID','Age','Sex','ErrorRate','Latency','AgeCategory'});
statsTable(:,"LunaID") = dataFirstVisit(:,"LunaID");
statsTable(:,"Age") = dataFirstVisit(:,"Age");
statsTable(:,"Sex") = dataFirstVisit(:,"Sex");
statsTable(:,"ErrorRate") = dataFirstVisit(:,"ErrorRate_noDropped_");
statsTable(:,"Latency") = dataFirstVisit(:,"AverageLatency_correctTrials_");
statsTable(:,"AgeCategory") = table(ageCategory);
writetable(statsTable,"FinalStatsData.xlsx")

%% effect size of post hoc analysis
meanC = 0.3054;
meanA = 0.2837;
meanAD = 0.2863;
varC = 0.00163; 
varA = 0.00177;
varAD = 0.00158;
stdC = 0.0427;
stdA = 0.0421;
stdAD = 0.0397;
nC = 34;
nA = 85;
nAD = 50;
pooledVarAC = ((nA-1)*varA + (nC-1)*varC)/(nA+nC-2);
pooledVarAAD = ((nA-1)*varA + (nAD-1)*varAD)/(nA+nAD-2);
pooledVarADC = ((nAD-1)*varAD + (nC-1)*varC)/(nAD+nC-2);
dfAC = nA+nC-2;
dfAAD = nA+nAD-2;
dfADC = nAD+nC-2;
tAC = (meanA-meanC)/(sqrt((pooledVarAC/nA)+(pooledVarAC/nC)));
tAAD = (meanA-meanAD)/(sqrt((pooledVarAAD/nA)+(pooledVarAAD/nAD)));
tADC = (meanAD-meanC)/(sqrt((pooledVarADC/nAD)+(pooledVarADC/nC)));
% effect size for pairwise comparison (correlation coefficient)
rAC = sqrt(tAC^2/(tAC^2 + dfAC));
rAAD = sqrt(tAAD^2/(tAAD^2 + dfAAD));
rADC = sqrt(tADC^2/(tADC^2 + dfADC));
% Cohen's d for pairwise comparison
dAC = abs(meanA-meanC)/pooledVarAC;
dAAD = abs(meanA-meanAD)/pooledVarAAD;
dADC = abs(meanAD-meanC)/pooledVarADC;
fprintf('Effect Size for Adults-Children: %1.4f\n',rAC)
fprintf('Effect Size for Adults-Adolescents: %1.4f\n',rAAD)
fprintf('Effect Size for Adolescents-Children: %1.4f\n',rADC)

%% repeated measures
clear
addpath('/Volumes/Hera/Abby/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
d = readtable('7t_eegAS_ErrorLatency.csv');
ids = unique(d.LunaID);
data = table('Size',[306 7],'VariableTypes',{'double','double','double','string','double','double','double'},...
    'VariableNames',{'LunaID','Scandate','Age','Sex','Latency','ErrorRate','VisitNumber'});
for i = 1:length(ids)
    ididx = find(d.LunaID==ids(i));
    for j = 1:length(ididx)
        subData = d(ididx(j),:);
        vistNumber = table(j);
        data(ididx(j),:) = [d(ididx(j),{'LunaID','ScanDate','Age','Sex','AverageLatency_correctTrials_','ErrorRate_noDropped_'}),vistNumber];
        %data(ididx(j),:) = [d(ididx(j),{'LunaID','Age','AverageLatency_correctTrials_','ErrorRate_noDropped_'}),vistNumber,sexID(subIDidx,"sex")];
    end
end
writetable(data,'RepeatedMeasuresLong.xlsx')

%% repeated measures anova
addpath('/Volumes/Hera/Abby/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/txt/')
d = readtable('7t_eegAS_ErrorLatency.csv');
ids = unique(d.LunaID);
data = table('Size',[306 7],'VariableTypes',{'double','double','double','string','double','double','double'},...
    'VariableNames',{'LunaID','ScanDate','Age','Sex','Latency','ErrorRate','VisitNumber'});

for i = 1:length(ids)
    ididx = find(d.LunaID==ids(i));
    for j = 1:length(ididx)
        subData = d(ididx(j),:);
        vistNumber = table(j);
        data(ididx(j),:) = [d(ididx(j),{'LunaID','ScanDate','Age','Sex','AverageLatency_correctTrials_','ErrorRate_noDropped_'}),vistNumber];
        %data(ididx(j),:) = [d(ididx(j),{'LunaID','Age','AverageLatency_correctTrials_','ErrorRate_noDropped_'}),vistNumber,sexID(subIDidx,"sex")];
    end
end

for k = 1:length(ids)
    subIdx = find(data.LunaID==ids(k));
    if length(subIdx) <3
        data(subIdx,:)=[];
    end
end
widedata = unstack(data,{'ScanDate','Age','Sex','Latency','ErrorRate'},'VisitNumber');
widedata(:,{'ErrorRate_x1','ErrorRate_x2','ErrorRate_x3','ScanDate_x1','ScanDate_x2','ScanDate_x3','Sex_x2','Sex_x3'})=[];
writetable(widedata,"RepeatedMeasuresWideFinalStats.xlsx")