# loading error latency table
errlatdata <- read.csv('/Volumes/Hera/Abby/AS_EEG/7t_eegAS_ErrorLatency_20250320.csv')
eegASdataCorrect <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/ID_info_correct_trials.csv')

# removing non-viable participants
library(dplyr)
library(LNCDR)
library(ggplot2)
viabledata <- inner_join(errlatdata,eegASdataCorrect,by = c('LunaID','ScanDate'))
colnames(viabledata)[colnames(viabledata)=='Age.x'] <- 'Age'

# histograms
hist(viabledata$Age)
hist(viabledata$AvgLatCor)
hist(viabledata$PercentErrCor)
hist(viabledata$CVLatCor)

# waterfall plot 
d <- data.frame(age=viabledata$Age, id=viabledata$LunaID,sex = viabledata$Sex)
p <- waterfall_plot(d)
plot(p)

# Average Correct Latency ~ 1/Age
invAgeLatModel <- lm(AvgLatCor ~ 1 + I(1/Age),data=viabledata)
summary(invAgeLatModel)
p <- ggplot(viabledata,aes(x=Age,y=AvgLatCor))+
  ggtitle('Average Correct Latency ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)

# CV Correct Latency ~ 1/Age
invAgeCVLatModel <- lm(CVLatCor ~ 1 + I(1/Age),data=viabledata)
summary(invAgeCVLatModel)
p <- ggplot(viabledata,aes(x=Age,y=CVLatCor))+
  ggtitle('Coefficient of Latency Variation (Correct Trials) ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)

# Percent Error Corrected ~ 1/Age
invAgePercentErrCorModel <- lm(PercentErrCor ~ 1 + I(1/Age),data=viabledata)
summary(invAgePercentErrCorModel)
p <- ggplot(viabledata,aes(x=Age,y=PercentErrCor))+
  ggtitle('Percent Error Corrected ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))
lunaize(p)

# Ignore stuff below this!!!
# Correct Trials Behavioral Data
# loading error rate and latency table
errorlatencydata <- read.csv('/Volumes/Hera/Abby/AS_EEG/7t_eegAS_ErrorLatency_20250310.csv')
eegASdataCorrect <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/CorrectTrials/ID_info_correct_trials.csv')

# removing non-viable participants
library(dplyr)
viabledata <- inner_join(errorlatencydata,eegASdataCorrect,by = c('LunaID','ScanDate'))
colnames(viabledata)[colnames(viabledata)=='Age.x'] <- 'Age'

# Looking at distributions of variables
hist(viabledata$Age)
hist(log1p(viabledata$ErrorRateBea))
hist(log1p(viabledata$ErrorRate))
hist(viabledata$AverageLatencyCorrect)
hist(viabledata$VarLatCor)

# Linear Models
library(ggplot2)
# lncd R package lunaize

# Latency
ageLatModel <- lm(AverageLatencyCorrect ~ 1 + Age, data=viabledata)
summary(ageLatModel)
AIC(ageLatModel)
ggplot(viabledata,aes(x=Age,y=AverageLatencyCorrect))+
  ggtitle('Average Correct Latency ~ 1 + Age')+
  geom_point() +
  geom_line(aes(group = LunaID),color = 'gray')+
  geom_smooth(method='lm',formula=y~1+x)

invAgeLatModel <- lm(AverageLatencyCorrect ~ 1 + I(1/Age), data=viabledata)
summary(invAgeLatModel)
AIC(invAgeLatModel)
ggplot(viabledata,aes(x=Age,y=AverageLatencyCorrect))+
  ggtitle('Average Correct Latency ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))

ageSexLatModel <- lm(AverageLatencyCorrect ~ 1 + Age + Sex, data=viabledata)
summary(ageSexLatModel)

invAgeSexLatModel <- lm(AverageLatencyCorrect ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexLatModel)

# Latency Variability
ageVarLatModel <- lm(VarLatCor ~ 1 + Age, data=viabledata)
summary(ageVarLatModel)
AIC(ageVarLatModel)
ggplot(viabledata,aes(x=Age,y=VarLatCor))+
  ggtitle('Correct Latency Variability ~ 1 + Age')+
  geom_point() +
  geom_line(aes(group = LunaID),color = 'gray')+
  geom_smooth(method='lm',formula=y~1+x)

invAgeVarLatModel <- lm(VarLatCor ~ 1 + I(1/Age), data=viabledata)
summary(invAgeVarLatModel)
AIC(invAgeVarLatModel)
ggplot(viabledata,aes(x=Age,y=VarLatCor))+
  ggtitle('Correct Latency Variability ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~I(1/x))

ageSexVarLatModel <- lm(VarLatCor ~ 1 + Age + Sex, data=viabledata)
summary(ageSexVarLatModel)

invAgeSexVarLatModel <- lm(VarLatCor ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexVarLatModel)

# Error Rate
ageErrRateModel <- lm(log1p(ErrorRate) ~ 1 + Age, data=viabledata)
summary(ageErrRateModel)
AIC(ageErrRateModel)
ggplot(viabledata,aes(x=Age,y=log1p(ErrorRate)))+
  ggtitle('Error Rate ~ Age Linear Regression')+
  geom_point()+
  geom_line(aes(group=LunaID),color='gray')+
  geom_smooth(method='lm',formula=log1p(y)~1+x)

invAgeErrRateModel <- lm(log1p(ErrorRate) ~ 1 + I(1/Age), data=viabledata)
summary(invAgeErrRateModel)
ggplot(viabledata,aes(x=1/Age,y=log1p(ErrorRate)))+
  ggtitle('Error Rate ~ Inverse Age Linear Regression')+
  geom_point()+
  geom_line(aes(group=LunaID),color='gray')+
  geom_smooth(method='lm',formula=log1p(y)~1+x)

ageSexErrRateModel <- lm(log1p(ErrorRate) ~ 1 + Age + Sex, data=viabledata)
summary(ageSexErrRateModel)

invAgeSexErrRateModel <- lm(log1p(ErrorRate) ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexErrRateModel)

# Incorrect Trials
# loading error rate and latency table
errorlatencydata <- read.csv('/Volumes/Hera/Abby/AS_EEG/7t_eegAS_ErrorLatency_20250310.csv')
eegASdataIncorrect <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/IncorrectTrials/ID_info_incorrect_trials.csv')

# removing non-viable participants
library(dplyr)
viabledata <- inner_join(errorlatencydata,eegASdataIncorrect,by = c('LunaID','ScanDate'))
colnames(viabledata)[colnames(viabledata)=='Age.x'] <- 'Age'

# Looking at distributions of variables
hist(viabledata$Age)
hist(viabledata$AvgerageLatencyIncorrect)
hist(viabledata$VarLatIncor)

# Linear Models
library(ggplot2)

# Latency
ageLatModel <- lm(AvgerageLatencyIncorrect ~ 1 + Age, data=viabledata)
summary(ageLatModel)
ggplot(viabledata,aes(x=Age,y=AvgerageLatencyIncorrect))+
  ggtitle('Average Incorrect Latency ~ 1 + Age')+
  geom_point() +
  geom_line(aes(group = LunaID),color = 'gray')+
  geom_smooth(method='lm',formula=y~1+x)

invAgeLatModel <- lm(AvgerageLatencyIncorrect ~ 1 + I(1/Age), data=viabledata)
summary(invAgeLatModel)
ggplot(viabledata,aes(x=1/Age,y=AvgerageLatencyIncorrect))+
  ggtitle('Average Incorrect Latency ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~1+x)

ageSexLatModel <- lm(AvgerageLatencyIncorrect ~ 1 + Age + Sex, data=viabledata)
summary(ageSexLatModel)

invAgeSexLatModel <- lm(AvgerageLatencyIncorrect ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexLatModel)

# Latency Variability
ageVarLatModel <- lm(VarLatIncor ~ 1 + Age, data=viabledata)
summary(ageVarLatModel)
ggplot(viabledata,aes(x=Age,y=VarLatIncor))+
  ggtitle('Incorrect Latency Variability ~ 1 + Age')+
  geom_point() +
  geom_line(aes(group = LunaID),color = 'gray')+
  geom_smooth(method='lm',formula=y~1+x)

invAgeVarLatModel <- lm(VarLatIncor ~ 1 + I(1/Age), data=viabledata)
summary(invAgeVarLatModel)
ggplot(viabledata,aes(x=1/Age,y=VarLatIncor))+
  ggtitle('Incorrect Latency Variability ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~1+x)

ageSexVarLatModel <- lm(VarLatIncor ~ 1 + Age + Sex, data=viabledata)
summary(ageSexVarLatModel)

invAgeSexVarLatModel <- lm(VarLatIncor ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexVarLatModel)

# Error Corrected Trials
# loading error rate and latency table
errorlatencydata <- read.csv('/Volumes/Hera/Abby/AS_EEG/7t_eegAS_ErrorLatency_20250310.csv')
eegASdataErrcorrect <- read.csv('/Volumes/Hera/Abby/AS_EEG/PrepPeriodAnalysis/ErrorCorrectTrials/ID_info_errorcorrected_trials.csv')

# removing non-viable participants
library(dplyr)
viabledata <- inner_join(errorlatencydata,eegASdataErrcorrect,by = c('LunaID','ScanDate'))
colnames(viabledata)[colnames(viabledata)=='Age.x'] <- 'Age'

# Looking at distributions of variables
hist(viabledata$Age)
hist(viabledata$AverageLatencyErrorCorrect)
hist(viabledata$VarLatErrCor)

# Linear Models
library(ggplot2)

# Latency
ageLatModel <- lm(AverageLatencyErrorCorrect ~ 1 + Age, data=viabledata)
summary(ageLatModel)
ggplot(viabledata,aes(x=Age,y=AverageLatencyErrorCorrect))+
  ggtitle('Average Error Corrected Latency ~ 1 + Age')+
  geom_point() +
  geom_line(aes(group = LunaID),color = 'gray')+
  geom_smooth(method='lm',formula=y~1+x)

invAgeLatModel <- lm(AverageLatencyErrorCorrect ~ 1 + I(1/Age), data=viabledata)
summary(invAgeLatModel)
ggplot(viabledata,aes(x=1/Age,y=AverageLatencyErrorCorrect))+
  ggtitle('Average Error Corrected Latency ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~1+x)

ageSexLatModel <- lm(AverageLatencyErrorCorrect ~ 1 + Age + Sex, data=viabledata)
summary(ageSexLatModel)

invAgeSexLatModel <- lm(AverageLatencyErrorCorrect ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexLatModel)

# Latency Variability
ageVarLatModel <- lm(VarLatErrCor ~ 1 + Age, data=viabledata)
summary(ageVarLatModel)
ggplot(viabledata,aes(x=Age,y=VarLatErrCor))+
  ggtitle('Error Corrected Latency Variability ~ 1 + Age')+
  geom_point() +
  geom_line(aes(group = LunaID),color = 'gray')+
  geom_smooth(method='lm',formula=y~1+x)

invAgeVarLatModel <- lm(VarLatErrCor ~ 1 + I(1/Age), data=viabledata)
summary(invAgeVarLatModel)
ggplot(viabledata,aes(x=1/Age,y=VarLatErrCor))+
  ggtitle('Error Corrected Latency Variability ~ 1 + 1/Age')+
  geom_point()+
  geom_line(aes(group = LunaID),color='gray')+
  geom_smooth(method='lm',formula=y~1+x)

ageSexVarLatModel <- lm(VarLatErrCor ~ 1 + Age + Sex, data=viabledata)
summary(ageSexVarLatModel)

invAgeSexVarLatModel <- lm(VarLatErrCor ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexVarLatModel)
