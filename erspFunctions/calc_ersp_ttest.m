function [b_intercept,t,p] = calc_ersp_ttest(T,erspdata,numTimes,numFreqs)
for i=1:numFreqs
    for j = 1:numTimes
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
end
