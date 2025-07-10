function [b,t] = calc_ersp_groupact(T,erspdata,numTimes,numFreqs)
% erspdata is (subjects,freqs,times)
% setting up waitbar
f = waitbar(0,'Starting','Name','Group Activation');
for i=1:numFreqs
    waitbar(i/numFreqs,f,sprintf('%d out of %d frequency points',i,numFreqs))
    for j = 1:numTimes
        % time,freq point data across subjects
        data = squeeze(erspdata(:,i,j));
        T.ersp = data;
        lme = fitlme(T,'ersp ~ 1 + (1 | id)');
        coeftable = lme.Coefficients;
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b(i,j) = double(coeftable(1,2));
        t(i,j) = double(coeftable(1,4));
    end
end
close(f)
end
