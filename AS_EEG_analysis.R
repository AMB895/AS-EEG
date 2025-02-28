# loading error rate and latency table
errorlatencydata <- read.csv('/Volumes/Hera/Abby/AS_EEG/7t_eegAS_ErrorLatency_01302025.csv')

# removing non-viable participants
viabledata <- errorlatencydata[errorlatencydata$Viable %in% 1,]

# Looking at distributions of variables
hist(viabledata$Age)
hist(log1p(viabledata$ErrorRateBea))
hist(log1p(viabledata$ErrorRate))
hist(viabledata$AverageLatencyCorrect)

# calculating accuracy 
viabledata['CorrAcc'] <- viabledata$Correct/(viabledata$Total-viabledata$Dropped)
viabledata['IncorrAcc'] <- viabledata$Incorrect/(viabledata$Total-viabledata$Dropped)
viabledata['ErrCorrAcc'] <- viabledata$ErrorCorrected/(viabledata$Total-viabledata$Dropped)
viabledata['TotalErrAcc'] <- (viabledata$ErrorCorrected+viabledata$Incorrect)/(viabledata$Total-viabledata$Dropped)

# Linear Models
library(ggplot2)

# Latency
ageLatModel <- lm(AverageLatencyCorrect ~ 1 + Age, data=viabledata)
summary(ageLatModel)
ggplot(viabledata,aes(x=Age,y=AverageLatencyCorrect))+
  geom_smooth(method='lm',formula=y~1+x)+
  ggtitle('Latency ~ Age Linear Regression')+
  geom_point()

invAgeLatModel <- lm(AverageLatencyCorrect ~ 1 + I(1/Age), data=viabledata)
summary(invAgeLatModel)
ggplot(viabledata,aes(x=1/Age,y=AverageLatencyCorrect))+
  geom_smooth(method='lm',formula=y~1+x)+
  ggtitle('Latency ~ Inverse Age Linear Regression')+
  geom_point()

ageSexLatModel <- lm(AverageLatencyCorrect ~ 1 + Age + Sex, data=viabledata)
summary(ageSexLatModel)

invAgeSexLatModel <- lm(AverageLatencyCorrect ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexLatModel)

# Error Rate
ageErrRateModel <- lm(log1p(ErrorRate) ~ 1 + Age, data=viabledata)
summary(ageErrRateModel)
ggplot(viabledata,aes(x=Age,y=log1p(ErrorRate)))+
  geom_smooth(method='lm',formula=log1p(y)~1+x)+
  ggtitle('Error Rate ~ Age Linear Regression')+
  geom_point()

invAgeErrRateModel <- lm(log1p(ErrorRate) ~ 1 + I(1/Age), data=viabledata)
summary(invAgeErrRateModel)
ggplot(viabledata,aes(x=1/Age,y=log1p(ErrorRate)))+
  geom_smooth(method='lm',formula=log1p(y)~1+x)+
  ggtitle('Error Rate ~ Inverse Age Linear Regression')+
  geom_point()

ageSexErrRateModel <- lm(log1p(ErrorRate) ~ 1 + Age + Sex, data=viabledata)
summary(ageSexErrRateModel)

invAgeSexErrRateModel <- lm(log1p(ErrorRate) ~ 1 + I(1/Age) + Sex, data=viabledata)
summary(invAgeSexErrRateModel)

# Accuracy
errCorAgeModel <- lm(ErrCorrAcc ~ 1 + Age,data=viabledata)
summary(errCorAgeModel)
ggplot(viabledata,aes(x=Age,y=ErrCorrAcc))+
  geom_smooth(method='lm',formula=y~1+x)+
  ggtitle('Error Corrected Rate ~ Age Linear Regression')+
  geom_point()

corAgeModel <- lm(CorrAcc ~ 1 + Age,data=viabledata)
summary(corAgeModel)
ggplot(viabledata,aes(x=Age,y=CorrAcc))+
  geom_smooth(method='lm',formula=y~1+x)+
  ggtitle('Correct Rate ~ Age Linear Regression')+
  geom_point()

incorAgeModel <- lm(IncorrAcc ~ 1 + Age, data=viabledata)
summary(incorAgeModel)
ggplot(viabledata,aes(x=Age,y=IncorrAcc))+
  geom_smooth(method='lm',formula=y~1+x)+
  ggtitle('Incorrect Rate ~ Age Linear Regression')+
  geom_point()
