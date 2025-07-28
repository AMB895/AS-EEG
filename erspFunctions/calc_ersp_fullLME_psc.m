function [b,t,AIC] = calc_ersp_fullLME_psc(T,erspdata1,erspdata2,numTimes,numFreqs)
% ersp data is (subjects,freqs,times)
for i = 1:numFreqs
    for j = 1:numTimes
        data1 = squeeze(erspdata1(:,i,j));
        data2 = squeeze(erspdata2(:,i,j));
        T.ersp = [data1;data2];
        % ensure id and trial type are categorical variables
        if ~isa(T.id,'categorical')
            T.id = categorical(T.id);
        end
        
        if ~isa(T.trialtype,'cateforical')
            T.trialtype = categorical(T.trialtype);
        end
        % run lmer
        lme = fitlme(T,'ersp ~ age*trialtype + (1 | id)');
        coeftable = lme.Coefficients;
        interceptcoeftable = coeftable(1,:);
        agecoeftable = coeftable(2,:);
        trialtypecoeftable = coeftable(3,:);
        interactioncoeftable = coeftable(4,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b_age(i,j) = double(agecoeftable(1,2));
        b_trialtype(i,j) = double(trialtypecoeftable(1,2));
        b_interaction(i,j) = double(interactioncoeftable(1,2));
        b_intercept(i,j) = double(interceptcoeftable(1,2));
        
        t_age(i,j) = double(agecoeftable(1,4));
        t_trialtype(i,j) = double(trialtypecoeftable(1,4));
        t_interaction(i,j) = double(interactioncoeftable(1,4));
        t_intercept(i,j) = double(interceptcoeftable(1,4));
        % get AIC
        AIC(i,j) = lme.ModelCriterion.AIC;
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
end