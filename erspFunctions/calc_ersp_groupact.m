function [b_intercept,t,p] = calc_ersp_groupact(T,erspdata,numTimes,numFreqs)
totalIterations = numTimes*numFreqs;
% setting up waitbar
f = waitbar(0,'Starting','Name','Group Activation');
for i=1:numFreqs
    for j = 1:numTimes
        % getting current iteration
        currentIteration = i*j;
        % updating waitbar
        waitbar(currentIteration/totalIterations,f,sprintf('%d out of %d iterations',currentIteration,totalIterations))
        
        data = squeeze(erspdata(:,i,j));
        T.power = data;
        lme = fitlme(T,'power ~ 1 + (1 | id)');
        coeftable = lme.Coefficients;
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b_intercept(i,j) = double(coeftable(1,2));
        t(i,j) = double(coeftable(1,4));
        p(i,j) = double(coeftable(1,6));
%         [~,p(i,j),~,stats] = ttest(data);
%         t(i,j) = stats.tstat;
    end
end
close(f)
end
