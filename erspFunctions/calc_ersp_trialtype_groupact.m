function [b,t] = calc_ersp_trialtype_groupact(T,erspdata1,erspdata2,numTimes,numFreqs)
% setting up waitbar
f = waitbar(0,'Starting','Name','Group Activation by Trial Type');
for i = 1:numFreqs
    waitbar(i/numFreqs,f,sprintf('%d out of %d frequency points',i,numFreqs))
    for j = 1:numTimes
        data1 = squeeze(erspdata1(:,i,j));
        data2 = squeeze(erspdata2(:,i,j));
        T.ersp = [data1;data2];
        lme = fitlme(T,'ersp ~ trialtype + (1 | id)');
        coeftable = lme.Coefficients;
        interceptcoeftable = coeftable(1,:);
        trialtypecoeftable = coeftable(2,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b_trialtype(i,j) = double(trialtypecoeftable(1,2));
        b_intercept(i,j) = double(interceptcoeftable(1,2));
        
        t_trialtype(i,j) = double(trialtypecoeftable(1,4));
        t_intercept(i,j) = double(interceptcoeftable(1,4));
    end
end
% saving coefficient values to single structure
b.trialtype = b_trialtype;
b.intercept = b_intercept;
% saving t values to single structure
t.trialtype = t_trialtype;
t.intercept = t_intercept;
close(f)
end