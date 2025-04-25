function [b_age,t,p] = calc_ersp_invageeffects(T,erspdata,numTimes,numFreqs)
% setting up waitbar
f = waitbar(0,'Starting','Name','Inverse Age Effects');
for i = 1:numFreqs
    for j = 1:numTimes
        % updating waitbar
        waitbar(i/numFreqs,f,sprintf('%d out of %d frequency points',i,numFreqs))
        
        data = squeeze(erspdata(:,i,j));
        T.power = data;
        lme = fitlme(T,'power ~ invage + (1 | id)');
        coeftable = lme.Coefficients;
        agecoeftable = coeftable(2,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        % Flip sign of age coefficient because input was inverse age
        b_age(i,j) = -1*double(agecoeftable(1,2));
        t(i,j) = -1*double(agecoeftable(1,4));
        p(i,j) = double(agecoeftable(1,6));
    end
end
close(f)
end