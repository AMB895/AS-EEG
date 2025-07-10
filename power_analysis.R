library(dplyr)
library(LNCDR)
library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyr)
clusterpowerdata <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/clusterpower_corvgsinteraction.csv');
signalchangedata <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/signalchange_corvgsinteraction.csv');
# cluster 1
# total power
data_long_totc1pow <- clusterpowerdata %>% 
  select(ID,VisitNum,Age,totcorC1Pow,totvgsC1Pow) %>%
  pivot_longer(cols = c(totcorC1Pow,totvgsC1Pow),names_to = "TrialType", values_to="C1Power") %>%
  mutate(
    TrialType = sub("tot(.*)C1Pow","\\1",TrialType),
    TrialType = as.factor(TrialType))
totC1PowerModel <- lmerTest::lmer(C1Power ~ 1 + InvAge + TrialType + InvAge:TrialType + (1 | ID), data = data_long_totc1pow %>% mutate(InvAge = 1/Age))
summary(totC1PowerModel)
p <- ggplot(data_long_totc1pow ,aes(x=Age,y=C1Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Total Power (dB)',color='TrialType')+
  ggtitle('Gamma Cluster (43.7 Hz)')
lunaize(p)

# relative power
data_long_relc1pow <- pivot_longer(clusterpowerdata,cols=c(relcorC1Pow,relvgsC1Pow),names_to='TrialType',values_to='C1Power')
p <- ggplot(data_long_relc1pow %>% mutate(power_z = scale(C1Power)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C1Power,color=TrialType))+
  #geom_point(alpha=0.5)+
  #geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Relative Power (dB)',color='TrialType')+
  ggtitle('Gamma Cluster (43.7 Hz)')
lunaize(p)

# signal change from baseline
data_long_c1sigchange <- pivot_longer(signalchangedata,cols=c(corC1SigChange,vgsC1SigChange),names_to='TrialType',values_to='C1SigChange')
p <- ggplot(data_long_c1sigchange %>% mutate(power_z = scale(C1SigChange)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C1SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Gamma Cluster (43.7 Hz)')
lunaize(p)

# cluster 2
# total power
data_long_totc2pow <- clusterpowerdata %>% 
  select(ID,VisitNum,Age,totcorC2Pow,totvgsC2Pow) %>%
  pivot_longer(cols = c(totcorC2Pow,totvgsC2Pow),names_to = "TrialType", values_to="C2Power") %>%
  mutate(
    TrialType = sub("tot(.*)C2Pow","\\1",TrialType),
    TrialType = as.factor(TrialType))
totC2PowerModel <- lmerTest::lmer(C2Power ~ 1 + InvAge + TrialType + InvAge:TrialType + (1 | ID), data = data_long_totc2pow %>% mutate(InvAge = 1/Age))
summary(totC2PowerModel)
p <- ggplot(data_long_totc2pow, aes(x=Age,y=C2Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',formula = y~I(1/x))+
  labs(x='Age',y='Total Power (dB)',color='TrialType')+
  ggtitle('High Beta Cluster (24.7-30.4 Hz)')
lunaize(p)

# relative power
data_long_relc2pow <- pivot_longer(clusterpowerdata,cols=c(relcorC2Pow,relvgsC2Pow),names_to='TrialType',values_to='C2Power')
p <- ggplot(data_long_relc2pow %>% mutate(power_z = scale(C2Power)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C2Power,color=TrialType))+
  #geom_point(alpha=0.5)+
  #geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Relative Power (dB)',color='TrialType')+
  ggtitle('High Beta Cluster (24.7-30.4 Hz)')
lunaize(p)

# signal change from baseline
data_long_c2sigchange <- pivot_longer(signalchangedata,cols=c(corC2SigChange,vgsC2SigChange),names_to='TrialType',values_to='C2SigChange')
p <- ggplot(data_long_c2sigchange %>% mutate(power_z = scale(C2SigChange)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C2SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('High Beta Cluster (24.7-30.4 Hz)')
lunaize(p)

# cluster 3
# total power
data_long_totc3pow <- clusterpowerdata %>% 
  select(ID,VisitNum,Age,totcorC3Pow,totvgsC3Pow) %>%
  pivot_longer(cols = c(totcorC3Pow,totvgsC3Pow),names_to = "TrialType", values_to="C3Power") %>%
  mutate(
    TrialType = sub("tot(.*)C3Pow","\\1",TrialType),
    TrialType = as.factor(TrialType))
totC3PowerModel <- lmerTest::lmer(C3Power ~ 1 + InvAge + TrialType + InvAge:TrialType + (1 | ID), data = data_long_totc3pow %>% mutate(InvAge = 1/Age))
summary(totC3PowerModel)
p <- ggplot(data_long_totc3pow,aes(x=Age,y=C3Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',formula=y~I(1/x),se=FALSE)+
  labs(x='Age (years)',y='Total Power (dB)',color='TrialType')+
  ggtitle('Alpha Cluster (7.4-12.0 Hz)')
lunaize(p)

# relative power
data_long_relc3pow <- pivot_longer(clusterpowerdata,cols=c(relcorC3Pow,relvgsC3Pow),names_to='TrialType',values_to='C3Power')
p <- ggplot(data_long_c3pow %>% mutate(power_z = scale(C3Power)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C3Power,color=TrialType))+
  #geom_point(alpha=0.5)+
  #geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age (years)',y='Relative Power (dB)',color='TrialType')+
  ggtitle('Alpha Cluster (7.4-12.0 Hz)')
lunaize(p)

# signal change from baseline
data_long_c3sigchange <- pivot_longer(signalchangedata,cols=c(corC3SigChange,vgsC3SigChange),names_to='TrialType',values_to='C3SigChange')
p <- ggplot(data_long_c3sigchange %>% mutate(power_z = scale(C3SigChange)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C3SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Alpha Cluster (7.4-12.0 Hz)')
lunaize(p)

# cluster 4
# total power
data_long_totc4pow <- clusterpowerdata %>% 
  select(ID,VisitNum,Age,totcorC4Pow,totvgsC4Pow) %>%
  pivot_longer(cols = c(totcorC4Pow,totvgsC4Pow),names_to = "TrialType", values_to="C4Power") %>%
  mutate(
    TrialType = sub("tot(.*)C4Pow","\\1",TrialType),
    TrialType = as.factor(TrialType))
totC4PowerModel <- lmerTest::lmer(C4Power ~ 1 + InvAge + TrialType + InvAge:TrialType + (1 | ID), data = data_long_totc4pow %>% mutate(InvAge = 1/Age))
summary(totC4PowerModel)
p <- ggplot(data_long_totc4pow ,aes(x=Age,y=C4Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',formula=y~I(1/x))+
  labs(x='Age',y='Total Power (dB)',color='TrialType')+
  ggtitle('Delta Cluster (3.3-4.1 Hz)')
lunaize(p)

# relative power
data_long_relc4pow <- pivot_longer(clusterpowerdata,cols=c(relcorC4Pow,relvgsC4Pow),names_to='TrialType',values_to='C4Power')
p <- ggplot(data_long_relc4pow  %>% mutate(power_z = scale(C4Power)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C4Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Relative Power (dB)',color='TrialType')+
  ggtitle('Delta Cluster (3.3-4.1 Hz)')
lunaize(p)

# signal change from baseline
data_long_c4sigchange <- pivot_longer(signalchangedata,cols=c(corC4SigChange,vgsC4SigChange),names_to='TrialType',values_to='C4SigChange')
p <- ggplot(data_long_c4sigchange  %>% mutate(power_z = scale(C4SigChange)[,1]) %>% filter(abs(power_z)<3),aes(x=Age,y=C4SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Delta Cluster (3.3-4.1 Hz)')
lunaize(p)
#################################
# Correct AS vs. Error AS Interaction Effect
clusterpowerdata <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/clusterpower_corerrcorinteraction.csv');
signalchangedata <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/signalchange_corerrcorinteraction.csv');
# cluster 1
data_long_c1pow <- pivot_longer(clusterpowerdata,cols=c(corC1Pow,errcorC1Pow),names_to='TrialType',values_to='C1Power')
p <- ggplot(data_long_c1pow,aes(x=Age,y=C1Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Power (dB)',color='TrialType')+
  ggtitle('Cluster 1 (4.2-10.0 Hz)')
lunaize(p)

data_long_c1sigchange <- pivot_longer(signalchangedata,cols=c(corC1SigChange,errcorC1SigChange),names_to='TrialType',values_to='C1SigChange')
p <- ggplot(data_long_c1sigchange,aes(x=Age,y=C1SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Cluster 1 (4.2-10.0 Hz)')
lunaize(p)

# cluster 2
data_long_c2pow <- pivot_longer(clusterpowerdata,cols=c(corC2Pow,errcorC2Pow),names_to='TrialType',values_to='C2Power')
p <- ggplot(data_long_c2pow,aes(x=Age,y=C2Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Power (dB)',color='TrialType')+
  ggtitle('Cluster 2 (21.2-32.3 Hz)')
lunaize(p)

data_long_c2sigchange <- pivot_longer(signalchangedata,cols=c(corC2SigChange,errcorC2SigChange),names_to='TrialType',values_to='C2SigChange')
p <- ggplot(data_long_c2sigchange,aes(x=Age,y=C2SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Cluster 2 (21.2-32.3 Hz)')
lunaize(p)

# cluster 3
data_long_c3pow <- pivot_longer(clusterpowerdata,cols=c(corC3Pow,errcorC3Pow),names_to='TrialType',values_to='C3Power')
p <- ggplot(data_long_c3pow,aes(x=Age,y=C3Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Power (dB)',color='TrialType')+
  ggtitle('Cluster 3 (53.9-55.6 Hz)')
lunaize(p)

data_long_c3sigchange <- pivot_longer(signalchangedata,cols=c(corC3SigChange,errcorC3SigChange),names_to='TrialType',values_to='C3SigChange')
p <- ggplot(data_long_c3sigchange,aes(x=Age,y=C3SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Cluster 3 (53.9-55.6 Hz)')
lunaize(p)

# cluster 4
data_long_c4pow <- pivot_longer(clusterpowerdata,cols=c(corC4Pow,errcorC4Pow),names_to='TrialType',values_to='C4Power')
p <- ggplot(data_long_c4pow,aes(x=Age,y=C4Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Power (dB)',color='TrialType')+
  ggtitle('Cluster 4 (53.9-57.3 Hz)')
lunaize(p)

data_long_c4sigchange <- pivot_longer(signalchangedata,cols=c(corC4SigChange,errcorC4SigChange),names_to='TrialType',values_to='C4SigChange')
p <- ggplot(data_long_c4sigchange,aes(x=Age,y=C4SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Cluster 4 (53.9-57.3 Hz)')
lunaize(p)

# cluster 5
data_long_c5pow <- pivot_longer(clusterpowerdata,cols=c(corC5Pow,errcorC5Pow),names_to='TrialType',values_to='C5Power')
p <- ggplot(data_long_c5pow,aes(x=Age,y=C5Power,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='Power (dB)',color='TrialType')+
  ggtitle('Cluster 5 (4.1-5.2 Hz)')
lunaize(p)

data_long_c5sigchange <- pivot_longer(signalchangedata,cols=c(corC5SigChange,errcorC5SigChange),names_to='TrialType',values_to='C5SigChange')
p <- ggplot(data_long_c5sigchange,aes(x=Age,y=C5SigChange,color=TrialType))+
  geom_point(alpha=0.5)+
  geom_line(alpha=0.25,aes(group=interaction(ID,TrialType)))+
  geom_smooth(method='lm',se=FALSE)+
  labs(x='Age',y='% Signal Change from Baseline',color='TrialType')+
  ggtitle('Cluster 5 (4.1-5.2 Hz)')
lunaize(p)
########################
data <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/interaction_avgpower.csv')
data_long_theta <- pivot_longer(data,cols=c(CorAS_ThetaPower,VGS_ThetaPower),names_to='TrialType',values_to='ThetaPower')
p <- ggplot(data_long_theta,aes(x=Age,y=ThetaPower,color=TrialType))+
  geom_point()+
  geom_line(aes(group=ID),color='gray')+
  geom_smooth(method='lm',se = FALSE)+
  labs(x='Age',y='Theta Power',color='Trial Type')
lunaize(p)

p <- ggplot(data_long_theta,aes(x=Age,y=ThetaPower,color=TrialType))+
  geom_smooth(method='lm',se = FALSE)+
  labs(x='Age',y='Theta Power',color='Trial Type')
lunaize(p)

data_long_highbeta <- pivot_longer(data,cols=c(CorAS_HighBetaPower,VGS_HighBetaPower),names_to='TrialType',values_to='HighBetaPower')
p <- ggplot(data_long_highbeta,aes(x=Age,y=HighBetaPower,color=TrialType))+
  geom_point()+
  geom_line(aes(group=ID),color='gray')+
  geom_smooth(method='lm',se = FALSE)+
  labs(x='Age',y='High Beta Power',color='Trial Type')
lunaize(p)

p <- ggplot(data_long_highbeta,aes(x=Age,y=HighBetaPower,color=TrialType))+
  geom_smooth(method='lm',se = FALSE)+
  labs(x='Age',y='High Beta Power',color='Trial Type')
lunaize(p)
#########################################################################
# Looking at Power ~ AgeGroup + TrialType + AgeGroup:TrialType
# Age as categorical variable

# Theta
theta_csvfiles <- list.files(path = "/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/AvgFreqTimeTables/theta",pattern = "*.csv",full.names=TRUE)
thetafile_info <- file.info(theta_csvfiles)
theta_csvfiles <- theta_csvfiles[order(thetafile_info$mtime,decreasing = FALSE)]
# Initialize a data frame to store all p-values
theta_pvals_df <- data.frame(File = character(), AgeGroup_p = numeric(), 
                          TrialType_p = numeric(), Interaction_p = numeric(), 
                          Intercept_p = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(theta_csvfiles)){
  thetadata_table <- read.csv(theta_csvfiles[i])
  print(paste("Loaded:",theta_csvfiles[i]))
  thetadata_table[4] <- as.factor(thetadata_table$AgeGroup)
  thetadata_table[6] <- as.factor(thetadata_table$TrialType)
  thetaModel <- lmerTest::lmer(Power ~ 1 + AgeGroup + TrialType + AgeGroup:TrialType + (1 | ID),data=thetadata_table %>% mutate(power_z = scale(Power)[ ,1]) %>% filter(abs(power_z)<3))
  thetacoef_table <- coef(summary(thetaModel))

    # Extract p-values
  theta_intercept_p <- thetacoef_table[1, "Pr(>|t|)"]  # p-value for the intercept
  theta_AgeGroup_p <- thetacoef_table["AgeGroup1", "Pr(>|t|)"]  # p-value for AgeGroup
  theta_trialtype_p <- thetacoef_table["TrialType2", "Pr(>|t|)"]  # p-value for TrialType
  theta_interaction_p <- thetacoef_table["AgeGroup1:TrialType2", "Pr(>|t|)"]  # p-value for interaction
  
  # get the file name from the full path
  theta_filename <- basename(theta_csvfiles[i])
  
  # store p-values in the data frame
  theta_pvals_df <- rbind(theta_pvals_df, data.frame(File = theta_filename, 
                                                  AgeGroup_p = theta_AgeGroup_p, 
                                                  TrialType_p = theta_trialtype_p, 
                                                  Interaction_p = theta_interaction_p, 
                                                  Intercept_p = theta_intercept_p))
}

# Apply FDR correction to all p-values
theta_pvals_df$AgeGroup_p_adj <- p.adjust(theta_pvals_df$AgeGroup_p, method = "fdr")
theta_pvals_df$TrialType_p_adj <- p.adjust(theta_pvals_df$TrialType_p, method = "fdr")
theta_pvals_df$Interaction_p_adj <- p.adjust(theta_pvals_df$Interaction_p, method = "fdr")
theta_pvals_df$Intercept_p_adj <- p.adjust(theta_pvals_df$Intercept_p, method = "fdr")

write.csv(theta_pvals_df,"/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/theta_pvals_catage.csv",row.names=FALSE)

# Low Beta
lowbeta_csvfiles <- list.files(path = "/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/AvgFreqTimeTables/lowbeta",pattern = "*.csv",full.names=TRUE)
lowbetafile_info <- file.info(lowbeta_csvfiles)
lowbeta_csvfiles <- lowbeta_csvfiles[order(lowbetafile_info$mtime,decreasing = FALSE)]
# Initialize a data frame to store all p-values
lowbeta_pvals_df <- data.frame(File = character(), AgeGroup_p = numeric(), 
                             TrialType_p = numeric(), Interaction_p = numeric(), 
                             Intercept_p = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(lowbeta_csvfiles)){
  lowbetadata_table <- read.csv(lowbeta_csvfiles[i])
  print(paste("Loaded:",lowbeta_csvfiles[i]))
  lowbetadata_table[4] <- as.factor(lowbetadata_table$AgeGroup)
  lowbetadata_table[6] <- as.factor(lowbetadata_table$TrialType)
  lowbetaModel <- lmerTest::lmer(Power ~ 1 + AgeGroup + TrialType + AgeGroup:TrialType + (1 | ID),data=lowbetadata_table %>% mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3))
  lowbetacoef_table <- coef(summary(lowbetaModel))
  
  # Extract p-values
  lowbeta_intercept_p <- lowbetacoef_table[1, "Pr(>|t|)"]  # p-value for the intercept
  lowbeta_AgeGroup_p <- lowbetacoef_table["AgeGroup1", "Pr(>|t|)"]  # p-value for AgeGroup
  lowbeta_trialtype_p <- lowbetacoef_table["TrialType2", "Pr(>|t|)"]  # p-value for TrialType
  lowbeta_interaction_p <- lowbetacoef_table["AgeGroup1:TrialType2", "Pr(>|t|)"]  # p-value for interaction
  
  # get the file name from the full path
  lowbeta_filename <- basename(lowbeta_csvfiles[i])
  
  # store p-values in the data frame
  lowbeta_pvals_df <- rbind(lowbeta_pvals_df, data.frame(File = lowbeta_filename, 
                                                     AgeGroup_p = lowbeta_AgeGroup_p, 
                                                     TrialType_p = lowbeta_trialtype_p, 
                                                     Interaction_p = lowbeta_interaction_p, 
                                                     Intercept_p = lowbeta_intercept_p))
}

# Apply FDR correction to all p-values
lowbeta_pvals_df$AgeGroup_p_adj <- p.adjust(lowbeta_pvals_df$AgeGroup_p, method = "fdr")
lowbeta_pvals_df$TrialType_p_adj <- p.adjust(lowbeta_pvals_df$TrialType_p, method = "fdr")
lowbeta_pvals_df$Interaction_p_adj <- p.adjust(lowbeta_pvals_df$Interaction_p, method = "fdr")
lowbeta_pvals_df$Intercept_p_adj <- p.adjust(lowbeta_pvals_df$Intercept_p, method = "fdr")

write.csv(lowbeta_pvals_df,"/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/lowbeta_pvals_catage.csv",row.names=FALSE)

# High Beta
highbeta_csvfiles <- list.files(path = "/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/AvgFreqTimeTables/highbeta",pattern = "*.csv",full.names=TRUE)
highbetafile_info <- file.info(highbeta_csvfiles)
highbeta_csvfiles <- highbeta_csvfiles[order(highbetafile_info$mtime,decreasing = FALSE)]
# Initialize a data frame to store all p-values
highbeta_pvals_df <- data.frame(File = character(), AgeGroup_p = numeric(), 
                               TrialType_p = numeric(), Interaction_p = numeric(), 
                               Intercept_p = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(highbeta_csvfiles)){
  highbetadata_table <- read.csv(highbeta_csvfiles[i])
  print(paste("Loaded:",highbeta_csvfiles[i]))
  highbetadata_table[4] <- as.factor(highbetadata_table$AgeGroup)
  highbetadata_table[6] <- as.factor(highbetadata_table$TrialType)
  highbetaModel <- lmerTest::lmer(Power ~ 1 + AgeGroup + TrialType + AgeGroup:TrialType + (1 | ID),data=highbetadata_table %>% mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3))
  highbetacoef_table <- coef(summary(highbetaModel))
  
  # Extract p-values
  highbeta_intercept_p <- highbetacoef_table[1, "Pr(>|t|)"]  # p-value for the intercept
  highbeta_AgeGroup_p <- highbetacoef_table["AgeGroup1", "Pr(>|t|)"]  # p-value for AgeGroup
  highbeta_trialtype_p <- highbetacoef_table["TrialType2", "Pr(>|t|)"]  # p-value for TrialType
  highbeta_interaction_p <- highbetacoef_table["AgeGroup1:TrialType2", "Pr(>|t|)"]  # p-value for interaction
  
  # get the file name from the full path
  highbeta_filename <- basename(highbeta_csvfiles[i])
  
  # store p-values in the data frame
  highbeta_pvals_df <- rbind(highbeta_pvals_df, data.frame(File = highbeta_filename, 
                                                         AgeGroup_p = highbeta_AgeGroup_p, 
                                                         TrialType_p = highbeta_trialtype_p, 
                                                         Interaction_p = highbeta_interaction_p, 
                                                         Intercept_p = highbeta_intercept_p))
}

# Apply FDR correction to all p-values
highbeta_pvals_df$AgeGroup_p_adj <- p.adjust(highbeta_pvals_df$AgeGroup_p, method = "fdr")
highbeta_pvals_df$TrialType_p_adj <- p.adjust(highbeta_pvals_df$TrialType_p, method = "fdr")
highbeta_pvals_df$Interaction_p_adj <- p.adjust(highbeta_pvals_df$Interaction_p, method = "fdr")
highbeta_pvals_df$Intercept_p_adj <- p.adjust(highbeta_pvals_df$Intercept_p, method = "fdr")

write.csv(highbeta_pvals_df,"/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/highbeta_pvals_catage.csv",row.names=FALSE)

#########################################################################
# Looking at Power ~ invAge + TrialType + invAge:TrialType
# Age as continuous variable

# Theta
theta_csvfiles <- list.files(path = "/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/AvgFreqTimeTables/theta",pattern = "*.csv",full.names=TRUE)
thetafile_info <- file.info(theta_csvfiles)
theta_csvfiles <- theta_csvfiles[order(thetafile_info$mtime,decreasing = FALSE)]
# Initialize a data frame to store all p-values
theta_pvals_df <- data.frame(File = character(), invAge_p = numeric(), 
                             TrialType_p = numeric(), Interaction_p = numeric(), 
                             Intercept_p = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(theta_csvfiles)){
  thetadata_table <- read.csv(theta_csvfiles[i])
  print(paste("Loaded:",theta_csvfiles[i]))
  thetadata_table[6] <- as.factor(thetadata_table$TrialType)
  thetaModel <- lmerTest::lmer(Power ~ 1 + invAge + TrialType + invAge:TrialType + (1 | ID),
                               data=thetadata_table %>% mutate(power_z = scale(Power)[ ,1]) %>% filter(abs(power_z)<3) %>% mutate(invAge=1/Age))
  thetacoef_table <- coef(summary(thetaModel))

    # Extract p-values
  theta_intercept_p <- thetacoef_table[1, "Pr(>|t|)"]  # p-value for the intercept
  theta_invAge_p <- thetacoef_table["invAge", "Pr(>|t|)"]  # p-value for AgeGroup
  theta_trialtype_p <- thetacoef_table["TrialType2", "Pr(>|t|)"]  # p-value for TrialType
  theta_interaction_p <- thetacoef_table["invAge:TrialType2", "Pr(>|t|)"]  # p-value for interaction
  
  # get the file name from the full path
  theta_filename <- basename(theta_csvfiles[i])
  
  # store p-values in the data frame
  theta_pvals_df <- rbind(theta_pvals_df, data.frame(File = theta_filename, 
                                                     invAge_p = theta_invAge_p, 
                                                     TrialType_p = theta_trialtype_p, 
                                                     Interaction_p = theta_interaction_p, 
                                                     Intercept_p = theta_intercept_p))
}

# Apply FDR correction to all p-values
theta_pvals_df$invAge_p_adj <- p.adjust(theta_pvals_df$invAge_p, method = "fdr")
theta_pvals_df$TrialType_p_adj <- p.adjust(theta_pvals_df$TrialType_p, method = "fdr")
theta_pvals_df$Interaction_p_adj <- p.adjust(theta_pvals_df$Interaction_p, method = "fdr")
theta_pvals_df$Intercept_p_adj <- p.adjust(theta_pvals_df$Intercept_p, method = "fdr")

write.csv(theta_pvals_df,"/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/theta_pvals_invage.csv",row.names=FALSE)

# Low Beta
lowbeta_csvfiles <- list.files(path = "/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/AvgFreqTimeTables/lowbeta",pattern = "*.csv",full.names=TRUE)
lowbetafile_info <- file.info(lowbeta_csvfiles)
lowbeta_csvfiles <- lowbeta_csvfiles[order(lowbetafile_info$mtime,decreasing = FALSE)]
# Initialize a data frame to store all p-values
lowbeta_pvals_df <- data.frame(File = character(), invAge_p = numeric(), 
                               TrialType_p = numeric(), Interaction_p = numeric(), 
                               Intercept_p = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(lowbeta_csvfiles)){
  lowbetadata_table <- read.csv(lowbeta_csvfiles[i])
  print(paste("Loaded:",lowbeta_csvfiles[i]))
  lowbetadata_table[6] <- as.factor(lowbetadata_table$TrialType)
  lowbetaModel <- lmerTest::lmer(Power ~ 1 + invAge + TrialType + invAge:TrialType + (1 | ID),
                                 data=lowbetadata_table %>% mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3) %>% mutate(invAge=1/Age))
  lowbetacoef_table <- coef(summary(lowbetaModel))
  
  # Extract p-values
  lowbeta_intercept_p <- lowbetacoef_table[1, "Pr(>|t|)"]  # p-value for the intercept
  lowbeta_invAge_p <- lowbetacoef_table["invAge", "Pr(>|t|)"]  # p-value for AgeGroup
  lowbeta_trialtype_p <- lowbetacoef_table["TrialType2", "Pr(>|t|)"]  # p-value for TrialType
  lowbeta_interaction_p <- lowbetacoef_table["invAge:TrialType2", "Pr(>|t|)"]  # p-value for interaction
  
  # get the file name from the full path
  lowbeta_filename <- basename(lowbeta_csvfiles[i])
  
  # store p-values in the data frame
  lowbeta_pvals_df <- rbind(lowbeta_pvals_df, data.frame(File = lowbeta_filename, 
                                                         invAge_p = lowbeta_invAge_p, 
                                                         TrialType_p = lowbeta_trialtype_p, 
                                                         Interaction_p = lowbeta_interaction_p, 
                                                         Intercept_p = lowbeta_intercept_p))
}

# Apply FDR correction to all p-values
lowbeta_pvals_df$invAge_p_adj <- p.adjust(lowbeta_pvals_df$invAge_p, method = "fdr")
lowbeta_pvals_df$TrialType_p_adj <- p.adjust(lowbeta_pvals_df$TrialType_p, method = "fdr")
lowbeta_pvals_df$Interaction_p_adj <- p.adjust(lowbeta_pvals_df$Interaction_p, method = "fdr")
lowbeta_pvals_df$Intercept_p_adj <- p.adjust(lowbeta_pvals_df$Intercept_p, method = "fdr")

write.csv(lowbeta_pvals_df,"/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/lowbeta_pvals_invage.csv",row.names=FALSE)

# High Beta
highbeta_csvfiles <- list.files(path = "/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/AvgFreqTimeTables/highbeta",pattern = "*.csv",full.names=TRUE)
highbetafile_info <- file.info(highbeta_csvfiles)
highbeta_csvfiles <- highbeta_csvfiles[order(highbetafile_info$mtime,decreasing = FALSE)]
# Initialize a data frame to store all p-values
highbeta_pvals_df <- data.frame(File = character(), invAge_p = numeric(), 
                                TrialType_p = numeric(), Interaction_p = numeric(), 
                                Intercept_p = numeric(), stringsAsFactors = FALSE)

for (i in seq_along(highbeta_csvfiles)){
  highbetadata_table <- read.csv(highbeta_csvfiles[i])
  print(paste("Loaded:",highbeta_csvfiles[i]))
  highbetadata_table[6] <- as.factor(highbetadata_table$TrialType)
  highbetaModel <- lmerTest::lmer(Power ~ 1 + invAge + TrialType + invAge:TrialType + (1 | ID),
                                  data=highbetadata_table %>% mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3) %>% mutate(invAge=1/Age))
  highbetacoef_table <- coef(summary(highbetaModel))
  
  # Extract p-values
  highbeta_intercept_p <- highbetacoef_table[1, "Pr(>|t|)"]  # p-value for the intercept
  highbeta_invAge_p <- highbetacoef_table["invAge", "Pr(>|t|)"]  # p-value for AgeGroup
  highbeta_trialtype_p <- highbetacoef_table["TrialType2", "Pr(>|t|)"]  # p-value for TrialType
  highbeta_interaction_p <- highbetacoef_table["invAge:TrialType2", "Pr(>|t|)"]  # p-value for interaction
  
  # get the file name from the full path
  highbeta_filename <- basename(highbeta_csvfiles[i])
  
  # store p-values in the data frame
  highbeta_pvals_df <- rbind(highbeta_pvals_df, data.frame(File = highbeta_filename, 
                                                           invAge_p = highbeta_invAge_p, 
                                                           TrialType_p = highbeta_trialtype_p, 
                                                           Interaction_p = highbeta_interaction_p, 
                                                           Intercept_p = highbeta_intercept_p))
}

# Apply FDR correction to all p-values
highbeta_pvals_df$invAge_p_adj <- p.adjust(highbeta_pvals_df$invAge_p, method = "fdr")
highbeta_pvals_df$TrialType_p_adj <- p.adjust(highbeta_pvals_df$TrialType_p, method = "fdr")
highbeta_pvals_df$Interaction_p_adj <- p.adjust(highbeta_pvals_df$Interaction_p, method = "fdr")
highbeta_pvals_df$Intercept_p_adj <- p.adjust(highbeta_pvals_df$Intercept_p, method = "fdr")

write.csv(highbeta_pvals_df,"/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/highbeta_pvals_invage.csv",row.names=FALSE)

###########################################################################
# Average frequency-time cluster models
thetatable <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/theta_avgClusterTable.csv')
thetatable[6] <- as.factor(thetatable$TrialType)
thetaModel <- lmerTest::lmer(Power ~ 1 + invAge + TrialType + invAge:TrialType + (1 | ID), data=thetatable %>% mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3) %>% mutate(invAge=1/Age))
summary(thetaModel)
AIC(thetaModel)

highbetatable <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/highbeta_avgClusterTable.csv')
highbetatable[6] <- as.factor(highbetatable$TrialType)
highbetaModel <- lmerTest::lmer(Power ~ 1 + invAge + TrialType + invAge:TrialType + (1 | ID), data=highbetatable %>%  mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3) %>% mutate(invAge=1/Age))
summary(highbetaModel)
AIC(highbetaModel)

thetatablesplit <- split(thetatable,thetatable$TrialType)
thetatablecor <- thetatablesplit$'1'

p<-ggplot(thetatablecor,aes(x=Age,y=Power))+
  ggtitle('Theta Power (Correct Trials) ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = ID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)

thetatablecor$AvgLatCor <- datatable$AvgLatCor
powerLatModel <- lmerTest::lmer(AvgLatCor ~ 1 + Power + (1 | ID),data=thetatablecor %>%  mutate(power_z = scale(Power)[,1]) %>% filter(abs(power_z)<3))
summary(powerLatModel)

##################################################################

library(dplyr)
library(LNCDR)
library(ggplot2)
library(lme4)
library(lmerTest)
datatable = read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorLatencyTable_subs_w_cor_and_errcor.csv')

# Behavioral stats
d <- data.frame(id = datatable$LunaID,age=datatable$Age)
p <- waterfall_plot(d)
lunaize(p)

summary(datatable)

hist(datatable$Age)

# Average Correct Latency
invAgeLatModel <- lmerTest::lmer(AvgLatCor ~ invAge + (1 | LunaID),data=datatable %>% mutate(invAge=1/Age))
summary(invAgeLatModel)
p <- ggplot(datatable,aes(x=Age,y=AvgLatCor))+
  ggtitle('Average Correct Latency ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)

# Percent Error Corrected
invAgePerErrCorModel <- lmerTest::lmer(PercentErrCor ~ invAge + (1 | LunaID),data=datatable %>% mutate(invAge=1/Age))
summary(invAgePerErrCorModel)
p <- ggplot(datatable,aes(x=Age,y=PercentErrCor))+
  ggtitle('Percent Error Corrected ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)

# Coefficient of Latency Variation (correct Trials)
invAgeCVLatModel <- lmerTest::lmer(CVLatCor ~ invAge+ (1 | LunaID),data=datatable %>% mutate(invAge=1/Age))
summary(invAgeCVLatModel)
p <- ggplot(datatable,aes(x=Age,y=CVLatCor))+
  ggtitle('Coefficient of Latency Variation (Correct Trials) ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)