function [p,t] = calc_ersp_paired_ttest(erspdata1,erspdata2,numTimes,numFreqs)
for i = 1:numFreqs
    for j = 1:numTimes
        data1 = squeeze(erspdata1(:,i,j));
        data2 = squeeze(erspdata2(:,i,j));
        [~,p(i,j),~,stats] = ttest(data1,data2);
        t(i,j) = stats.tstat;
    end
end
end