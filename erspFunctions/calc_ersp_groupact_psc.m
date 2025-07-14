function [b,t] = calc_ersp_groupact_psc(T,erspdata,numTimes,numFreqs)
% erspdata is (subjects,freqs,times)
for i=1:numFreqs
    for j=1:numTimes
        % linear mixed effects model at each time frequency point
        data = squeeze(erspdata(:,i,j));
        T.ersp = data;
        lme = fitlme(T,'ersp ~ 1 + (1|id)');
        coeftable = lme.Coefficients;
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b(i,j) = double(coeftable(1,2));
        t(i,j) = double(coeftable(1,4));
    end
end
end