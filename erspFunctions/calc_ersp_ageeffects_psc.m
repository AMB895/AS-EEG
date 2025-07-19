function [b,t] = calc_ersp_ageeffects_psc(T,erspdata,numTimes,numFreqs)
% erspdata is (subjects,freqs,times)
for i = 1:numFreqs
    for j = 1:numTimes
         % time,freq point data across subjects
        data = squeeze(erspdata(:,i,j));
        T.ersp = data;
        lme = fitlme(T,'ersp ~ 1 + age + (1 | id)');
        coeftable = lme.Coefficients;
        agecoeftable = coeftable(2,:);
        intercepttable = coeftable(1,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b_age(i,j) = double(agecoeftable(1,2));
        b_intercept(i,j) = double(intercepttable(1,2));
        t_age(i,j) = double(agecoeftable(1,4));
        t_intercept(i,j) = double(intercepttable(1,4));
    end
end
% add b values to struct
b.age = b_age;
b.intercept = b_intercept;
% add t values to struct 
t.age = t_age;
t.intercept = t_intercept;
end