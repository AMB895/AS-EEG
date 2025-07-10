function [b,t] = calc_ersp_invageeffects(T,erspdata,numTimes,numFreqs)
% erspdata is (subjects,freqs,times)
% setting up waitbar
f = waitbar(0,'Starting','Name','Inverse Age Effects');
for i = 1:numFreqs
        % updating waitbar
        waitbar(i/numFreqs,f,sprintf('%d out of %d frequency points',i,numFreqs))
    for j = 1:numTimes
         % time,freq point data across subjects
        data = squeeze(erspdata(:,i,j));
        T.ersp = data;
        lme = fitlme(T,'ersp ~ 1 + invage + (1 | id)');
        coeftable = lme.Coefficients;
        agecoeftable = coeftable(2,:);
        intercepttable = coeftable(1,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        % Flip sign of age coefficient because input was inverse age
        b_age(i,j) = -1*double(agecoeftable(1,2));
        b_intercept(i,j) = double(intercepttable(1,2));
        t_age(i,j) = -1*double(agecoeftable(1,4));
        t_intercept(i,j) = double(intercepttable(1,4));
    end
end
% add b values to struct
b.age = b_age;
b.intercept = b_intercept;
% add t values to struct 
t.age = t_age;
t.intercept = t_intercept;
close(f)
end