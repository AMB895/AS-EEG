function [b,t] = calc_ersp_trialtype_groupact_psc(T,erspdata1,erspdata2,numTimes,numFreqs,ctrl4age)
% erspdata is (subjects,freqs,times)
for i = 1:numFreqs
    for j = 1:numTimes
        data1 = squeeze(erspdata1(:,i,j));
        data2 = squeeze(erspdata2(:,i,j));
        T.ersp = [data1;data2];
        if ctrl4age
            % control for linear age (not looking at developmental age effects)
            lme = fitlme(T,'ersp ~ trialtype + eeg_age + (1 | id)');
            coeftable = lme.Coefficients;
            interceptcoeftable = coeftable(1,:);
            agecoeftable = coeftable(2,:);
            trialtypecoeftable = coeftable(3,:);
            % Coefficient Tables with rows: 
    % Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
            % Flip sign of age coefficient because input was inverse age
            b_trialtype(i,j) = double(trialtypecoeftable(1,2));
            b_age(i,j) = -1*double(agecoeftable(1,2));
            b_intercept(i,j) = double(interceptcoeftable(1,2));

            t_trialtype(i,j) = double(trialtypecoeftable(1,4));
            t_age(i,j) = -1*double(agecoeftable(1,4));
            t_intercept(i,j) = double(interceptcoeftable(1,4));
        else
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
end

% saving coef values and t-values to single structure
if ctrl4age
    b.trialtype = b_trialtype;
    b.intercept = b_intercept;
    b.age = b_age;
    
    t.trialtype = t_trialtype;
    t.intercept = t_intercept;
    t.age = t_age;
else
    b.trialtype = b_trialtype;
    b.intercept = b_intercept;
    
    t.trialtype = t_trialtype;
    t.intercept = t_intercept;
end
end