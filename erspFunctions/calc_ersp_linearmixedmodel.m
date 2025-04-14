function [b,t,p] = calc_ersp_linearmixedmodel(T,erspdata1,erspdata2,numTimes,numFreqs)
totalIterations = numTimes*numFreqs;
% setting up waitbar
f = waitbar(0,'Starting','Name','Group Activation');
for i = 1:numFreqs
    for j = 1:numTimes
        % getting current iteration
        currentIteration = i*j;
        % updating waitbar
        waitbar(currentIteration/totalIterations,f,sprintf('%d out of %d iterations',currentIteration,totalIterations))
        
        data1 = squeeze(erspdata1(:,i,j));
        data2 = squeeze(erspdata2(:,i,j));
        T.power = [data1;data2];
        lme = fitlme(T,'power ~ invage + (1 | id)');
        coeftable = lme.Coefficients;
        agecoeftable = coeftable(2,:);
        trialtypecoeftable = coeftable(3,:);
        interactioncoeftable = coeftable(4,:);
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        % Flip sign of age coefficient because input was inverse age
        b_age(i,j) = -1*double(agecoeftable(1,2));
        b_trialtype(i,j) = double(trialtypecoeftable(1,2));
        b_interaction(i,j) = double(interactioncoeftable(1,2));
        t_age(i,j) = double(agecoeftable(1,4));
        t_trialtype(i,j) = double(trialtypecoeftable(1,4));
        t_interaction(i,j) = double(interactioncoeftable(1,4));
        p_age(i,j) = double(agecoeftable(1,6));
        p_trialtype(i,j) = double(trialtypecoeftable(1,6));
        p_interaction(i,j) = double(interactioncoeftable(1,6));
    end
end
% saving coefficient values to single structure
b.age = b_age;
b.trialtype = b_trialtype;
b.interaction = b_interaction;
% saving t values to single structure
t.age = t_age;
t.trialtype = t_trialtype;
t.interaction = t_interaction;
% saving p values to single structure
p.age = p_age;
p.trialtype = p_trialtype;
p.interaction = p_interaction;
close(f)
end