function [b,t,AIC] = calc_ersp_groupact_psc(T,erspdata,numTimes,numFreqs)
% erspdata is (subjects,freqs,times)
for i=1:numFreqs
    for j=1:numTimes
        % linear mixed effects model at each time frequency point
        data = squeeze(erspdata(:,i,j));
        T.ersp = data;
        % ensure id is categorical variable
        if ~isa(T.id,'categorical')
            T.id = categorical(T.id);
        end
        lme = fitlme(T,'ersp ~ 1 + (1|id)');
        coeftable = lme.Coefficients;
        % Coefficient Tables with rows: 
% Name, Estimate, SE, tStat, DF, pValue, Lower, Upper
        b(i,j) = double(coeftable(1,2));
        t(i,j) = double(coeftable(1,4));
        
        % Get AIC 
        AIC(i,j) = lme.ModelCriterion.AIC;
    end
end
end