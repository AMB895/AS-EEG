function [b_intercept,t,p] = calc_ersp_groupact(T,erspdata,numTimes,numFreqs)
% setting up waitbar
f = waitbar(0,'Starting','Name','Group Activation');
for i=1:numFreqs
    waitbar(i/numFreqs,f,sprintf('%d out of %d frequency points',i,numFreqs))
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
    end
end
close(f)
end
