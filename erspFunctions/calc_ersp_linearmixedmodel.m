function [b,t] = calc_ersp_linearmixedmodel(T,erspdata1,erspdata2,numTimes,numFreqs)
% setting up waitbar
f = waitbar(0,'Starting','Name','Linear Mixed Model');
for i = 1:numFreqs
    waitbar(i/numFreqs,f,sprintf('%d out of %d frequency points',i,numFreqs))
    for j = 1:numTimes
        data1 = squeeze(erspdata1(:,i,j));
        data2 = squeeze(erspdata2(:,i,j));
        T.ersp = [data1;data2];
        lme = fitlme(T,'ersp ~ invage*trialtype + (1 | id)');
        % can include people that have one or the other and not do weird
        % matching thing
        coeftable = lme.Coefficients;
        interceptcoeftable = coeftable(1,:);
        agecoeftable = coeftable(2,:);
        trialtypecoeftable = coeftable(3,:);
        interactioncoeftable = coeftable(4,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        % Flip sign of age coefficient because input was inverse age
        b_age(i,j) = -1*double(agecoeftable(1,2));
        b_trialtype(i,j) = double(trialtypecoeftable(1,2));
        b_interaction(i,j) = double(interactioncoeftable(1,2));
        b_intercept(i,j) = double(interceptcoeftable(1,2));
        
        % Flip sign of t stat because input was inverse age
        t_age(i,j) = -1*double(agecoeftable(1,4));
        t_trialtype(i,j) = double(trialtypecoeftable(1,4));
        t_interaction(i,j) = double(interactioncoeftable(1,4));
        t_intercept(i,j) = double(interceptcoeftable(1,4));
    end
end
% saving coefficient values to single structure
b.age = b_age;
b.trialtype = b_trialtype;
b.interaction = b_interaction;
b.intercept = b_intercept;
% saving t values to single structure
t.age = t_age;
t.trialtype = t_trialtype;
t.interaction = t_interaction;
t.intercept = t_intercept;
close(f)
end