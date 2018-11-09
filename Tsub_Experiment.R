library(chron)
library(dplyr)
library(ggplot2)
library(lme4)
source("HighStatLibV10.R")
library(lattice)
library(multcomp)
unloadNamespace("lmerTest")
library(lmerTest)
library(ReporteRs)
library(broom)

###########################################################################################Function for treatment validation

Valid<-function(modelobject){
Resid<- resid(modelobject,type='pearson')
Fitted<- fitted(modelobject)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)
qqnorm(Resid)
qqline(Resid)
}

###########################################################################################Effect of treatment on T values and Treatment validation

#####Subset the dataframe so that there is only the info that I need for the paper
Plasma_T <-read.csv("Tsub_PlasmaT_20160702.csv") %>%  
  select(SampleID,AnimalID,TreatmentAll,Colony,Testosterone) %>% 
  filter(TreatmentAll!="GnRH",TreatmentAll!="Post-Conflict") %>% 
  rename(Treatment=TreatmentAll)
View(Plasma_T)
levels(Plasma_T$Treatment)

write.csv(Plasma_T,"TreatmentEffect_Plasma.csv",row.names = FALSE)


#####model: Treatment and breeding status effect on plasma T must change name to breeding female
T_Plasma<-read.csv("TreatmentEffect_Plasma.csv") %>%   
  mutate(Treatment=ifelse(Treatment=="Queen","Breeding female",as.character(Treatment))) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Control", "Tlow","Thigh","Breeding female"))) 
levels(T_Plasma$Treatment)
View(T_Plasma)

T_Effect1<-lmer(log(Testosterone)~Treatment+(1|AnimalID)+(1|Colony),data=T_Plasma)
print(summary(T_Effect1),digits=3)


T_Effect2<-update(T_Effect1,~.-Treatment)
summary(T_Effect2)
anova(T_Effect1,T_Effect2)
drop1(T_Effect1,test="Chisq")


####post hoc comparison: computed all possible comparison
T_Plasma_ph <- glht(T_Effect1, linfct = mcp(Treatment='Tukey'))
print(summary(T_Plasma_ph),digits=3)

####Graph
Fig1<-ggplot(T_Plasma,aes(Treatment,Testosterone))+
  geom_boxplot(aes(fill=Treatment),outlier.shape=NA)+
  geom_point(size=3,alpha=0.4)+
  #geom_jitter(width=0,height=0.2,size=2,col=1,alpha=0.6)+
  geom_line(aes(group=AnimalID),alpha=0.3)+
  scale_x_discrete("",expand=c(0.05,0.05))+
  theme(axis.text.x = element_text(size = 12))+
  scale_y_continuous("Plasma Testosterone [ng/dl]",limits=c(0,0.9))+
  theme_classic()+
  scale_fill_manual(values = c("white", "#00CCFF","#FF3300", "red"))+
  guides(fill=FALSE)



####Model Validation

Resid<- resid(T_Effect1,type="pearson")
Fitted<- fitted(T_Effect1)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)
#The gaussian model shows a lot of heterogeneity, Also there are 4 negative fitted values
#It looks much better with the log transformed data 
#the gamma model is good as well 


#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)
#Gaussian: Not bad at all 


###########################################################################################Effect of treatment on T values and Treatment on urine levels 


####Data preparation must change name to breeding female
T_Urine<-read.csv("T_UrineQueen.csv")%>% 
  mutate(TreatmentAll=ifelse(TreatmentAll=="Queen","Breeding female",as.character(TreatmentAll))) %>% 
  mutate(TreatmentAll = factor(TreatmentAll, levels = c("No T", "Tlow","Thigh","Breeding female"))) %>% 
  mutate(Treatment = factor(Treatment, levels = c("Baseline", "Control", "Tlow","Thigh","Queen"))) %>% 
  filter(AnimalID!='G3F008') %>% 
  na.omit()
names(T_Urine)
levels(T_Urine$Treatment)
levels(T_Urine$TreatmentAll)
table (T_Urine$Treatment)
  
####Model where pooled baseline and control 

T_Effect_U1<-lmer(log(Testosterone)~TreatmentAll+(1|DTB_Colony/AnimalID),data=T_Urine)
drop1(T_Effect_U1,test="Chisq")
summary(T_Effect_U1,digits=3)

#write model output
tidy(T_Effect_U1) %>% 
  write.csv(.,'Tvalidation_Urine.csv')


T_Urine_ph <- glht(T_Effect_U1, linfct = mcp(TreatmentAll='Tukey'))
T_Urine_ph_Summary<-print(summary(T_Urine_ph),digits=3)
tapply(T_Urine$Weight,T_Urine$Treatment,mean)


####Model Validation
Resid<- resid(T_Effect_U1,type="pearson")
Fitted<- fitted(T_Effect_U1)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)
#The gaussian model shows a lot of heterogeneity, Also there are 4 negative fitted values
#It looks much better with the log transformed data 
#the gamma model is good as well



#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)
#Gaussian: Not that great, not that bad

####Graph
scaleFUN <- function(x) sprintf("%.1f", x)


Fig2<-ggplot(T_Urine,aes(TreatmentAll,Testosterone))+
  geom_boxplot(aes(fill=TreatmentAll),outlier.shape=NA)+
  geom_point(size=3,alpha=0.3)+
  #geom_jitter(width=0,height=0.2,size=2,col=1,alpha=0.6)+
  geom_line(aes(group=AnimalID),alpha=0.3)+
  scale_x_discrete("",expand=c(0.05,0.05))+
  theme(axis.text.x = element_text(size = 12))+
  scale_y_continuous("Urine Testosterone [ng/dl]",limits=c(0,70),labels=scaleFUN)+
  theme_classic()+
  scale_fill_manual(values = c("grey", "#00CCFF","#FF3300", "red"))+
  guides(fill=FALSE)

###########################################################################################Effect of treatment on helping

Focal_Treatment<-read.csv("Focal_Treatment_20160908.csv") %>% 
   mutate(Treatment = factor(Treatment, levels = c( "Control", "Tlow","Thigh"))) %>%
  mutate(CenteredDur=as.numeric(scale(DurationSeen,scale=FALSE)),ScaledDur=as.numeric(scale(DurationSeen)),Burrow=Dig+Sweep+Kick,HelpRate=WorkNoGnaw/DurationSeen,WorkRate=WorkNoGnaw/DurationSeen,BurrowRate=Burrow/DurationSeen,FCRate=FoodCarry/DurationSeen,NRate=Nest/DurationSeen) 



############log transformed Treatment:Nothing. Generally some of the gamma model are not better, some don t converge and give unrealistic p-value
Work1_T<-lmer(log(WorkNoGnaw)~Treatment+(1|Colony/AnimalID),data=Focal_Treatment)
summary(Work1_T)
tidy(Work1_T)
drop1(Work1_T,test="Chisq")

Work1_T<-glmer(WorkNoGnaw~Treatment+(1|Colony/AnimalID),data=Focal_Treatment,family="gaussian"(link='log'))
summary(Work1_T)
tidy(Work1_T)
drop1(Work1_T,test="Chisq")

Work1_T<-lmer(WorkNoGnaw~Treatment+(1|Colony/AnimalID),data=Focal_Treatment)
summary(Work1_T)
tidy(Work1_T)
drop1(Work1_T,test="Chisq")

tidy(Work1_T) %>% 
  write.csv(.,'T_help.csv')

glht(Work1_T, linfct = mcp(Treatment='Tukey'))



Resid<- resid(Work1_T,type="pearson")
Fitted<- fitted(Work1_T)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

Fig3<-ggplot(Focal_Treatment,aes(Treatment,WorkNoGnaw))+
  geom_boxplot(aes(fill=Treatment),outlier.shape=NA)+
  geom_point(size=3,alpha=0.3)+
  #geom_jitter(width=0,height=0.2,size=2,col=1,alpha=0.6)+
  geom_line(aes(group=AnimalID),alpha=0.3)+
  scale_x_discrete("",expand=c(0.05,0.05))+
  theme(axis.text.x = element_text(size = 12))+
  scale_y_continuous("Hepling duration in min",limits=c(0,120))+
  theme_classic()+
  scale_fill_manual(values = c("white", "#00CCFF","#FF3300"))+
  guides(fill=FALSE)



###############################################################################################Effect on Burrow
Burrow_T<-lmer(log(Burrow)~Treatment+(1|Colony/AnimalID),data=Focal_Treatment)
summary(Burrow_T)
tidy(Burrow_T)
drop1(Burrow_T,test="Chisq")

tidy(Burrow_T) %>% 
  write.csv(.,'T_Burrow.csv')

summary(glht(Burrow_T, linfct = mcp(Treatment='Tukey')))

Burrow_T<-glmer(Burrow~Treatment+(1|Colony/AnimalID),family="gaussian"(link='log'),data=Focal_Treatment)
summary(Burrow_T)




Resid<- resid(Burrow_T,type="pearson")
Fitted<- fitted(Burrow_T)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)


###############################################################################################Effect on Nest Building 

Nest_T<-lmer(log(Nest)~Treatment+(1|Colony/AnimalID),data=Focal_Treatment %>% 
    filter(Nest>0))
summary(Nest_T)
tidy(Nest_T)
drop1(Nest_T,test="Chisq")

tidy(Nest_T) %>% 
  write.csv(.,'T_Nest.csv')

summary(glht(Nest_T, linfct = mcp(Treatment='Tukey')))



Resid<- resid(Nest_T,type="pearson")
Fitted<- fitted(Nest_T)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)


###############################################################################################Effect on FoodCarry 

FC_T<-lmer(log(FoodCarry)~Treatment+(1|Colony/AnimalID),data=Focal_Treatment %>% 
    filter(FoodCarry>0))
summary(FC_T)
tidy(FC_T)
drop1(FC_T,test="Chisq")

tidy(FC_T) %>% 
  write.csv(.,'T_FC.csv')

summary(glht(FC_T, linfct = mcp(Treatment='Tukey')))



Resid<- resid(FC_T,type="pearson")
Fitted<- fitted(FC_T)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)




###############################################################################################Relative Helping, I must work with the scale duration. The response variable can be left intact, log transformed, or a gamma model can be used. I have decided to use log transformed data
 
######### non transformed data, non homogeneity of variace, no quadratic effect 
Work_TR<-lmer(WorkNoGnaw~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Work_TR),digits=3)
drop1(Work_TR,test='Chisq')

######### ln transformed data
Work_TR<-lmer(log(WorkNoGnaw)~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Work_TR),digits=3)
drop1(Work_TR,test='Chisq')

tidy(Work_TR) %>% 
  write.csv(.,'Work_Relative.csv')


Work_TR <- glht(Work_TR, linfct = mcp(Treatment='Tukey'))
print(summary(Work_TR),digits=3)



######### log link, pattern in residual 
Work_TR<-glmer(WorkNoGnaw~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID), family="gaussian"(link='log'),data=Focal_Treatment)
print(summary(Work_TR),digits=3)
drop1(Work_TR,test='Chisq')


######### Gamma. there is a pattern in the residual, 
Work_TR<-glmer(WorkNoGnaw~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=Gamma(link=log),data=Focal_Treatment)
print(summary(Work_TR),digits=3)
drop1(Work_TR,test='Chisq')


######### Rate with beta family
poop<-glmmTMB(WorkRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=list(family="beta",link="logit"),data=Focal_Treatment)
summary(Work_TR)
resid(Work_TR)



##########Data visualisation about the relationship

ggplot(Focal_Treatment,aes(CenteredDur,HelpRate))+
  geom_point()+
  geom_smooth(method='auto')

ggplot(Focal_Treatment,aes(ScaledDur,(Kick+Dig+Sweep)))+
  geom_point()+
  geom_smooth(method='auto')



par(mfrow=c(1,1))

Resid<- resid(Work_TR,type='pearson')
Fitted<- fitted(Work_TR)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

###############################################################################################Relative Burrowing 
 
######### non transformed data, non homogeneity of variace, no quadratic effect 
Burrow_TR<-lmer(Burrow~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Burrow_TR),digits=3)
drop1(Burrow_TR,test='Chisq')

######### ln transformed data, almost opposite funnel shape smaller variance for larger fitted value 
Burrow_TR<-lmer(log(Burrow)~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Burrow_TR),digits=3)
drop1(Burrow_TR,test="Chisq")


tidy(Burrow_TR) %>% 
  write.csv(.,'T_Burrow_Relative.csv')


Burrow_TR_PH <- glht(Burrow_TR, linfct = mcp(Treatment='Tukey'))
summary(Burrow_TR_PH)


######### log link, not ideal
Burrow_TR<-glmer(Burrow~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID), family="gaussian"(link='log'),data=Focal_Treatment)
print(summary(Burrow_TR),digits=3)
drop1(Burrow_TR,test='Chisq')


######### Gamma. there is a pattern in the residual, 
Burrow_TR<-glmer(WorkNoGnaw~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=Gamma(link=log),data=Focal_Treatment)
print(summary(Burrow_TR),digits=3)
drop1(Burrow_TR,test='Chisq')


######### Rate 
Burrow_TR<-lmer(BurrowRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Burrow_TR),digits=3)
drop1(Burrow_TR,test='Chisq')


######### Rate with beta family, distribution of residuals is not that awesome
Burrow_TR<-glmmTMB(BurrowRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=list(family="beta",link="logit"),data=Focal_Treatment)


##########Data visualisation about the relationship

ggplot(Focal_Treatment,aes(CenteredDur,HelpRate))+
  geom_point()+
  geom_smooth(method='auto')

ggplot(Focal_Treatment,aes(ScaledDur,(Kick+Dig+Sweep)))+
  geom_point()+
  geom_smooth(method='auto')



par(mfrow=c(1,1))

Resid<- resid(Burrow_TR,type='pearson')
Fitted<- fitted(Burrow_TR)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

Resid<- resid(Burrow_TR)
Fitted<- fitted(Burrow_TR)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)


###############################################################################################Relative Food Nest Building 


######### non transformed data, non homogeneity of variace, no quadratic effect 
Nest_TR<-lmer(Nest~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')

table(Focal_Treatment %>% filter(Nest>0) %>% select(Treatment))

######### ln transformed data, all good 
Nest_TR<-lmer(log(Nest)~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment %>% filter(Nest>0))
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')

tidy(Nest_TR) %>% 
  write.csv(.,'T_Nest_FC_Relative.csv')


Nest_TR_PH <- glht(Nest_TR, linfct = mcp(Treatment='Tukey'))
summary(Nest_TR_PH)


######### log link, not ideal not super homogenous 
Nest_TR<-glmer(Nest~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID), family="gaussian"(link='log'),data=Focal_Treatment %>% filter(Nest>0))
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')


######### Gamma. quite homogenous  
Nest_TR<-glmer(Nest~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=Gamma(link=log),data=Focal_Treatment %>% filter(Nest>0))
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')


######### Rate: terrible
Nest_TR<-lmer(NRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')


######### Rate with beta family
poop<-glmmTMB(NRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=list(family="beta",link="logit"),data=Focal_Treatment %>% filter(Nest>0))

par(mfrow=c(1,1))

Resid<- resid(Nest_TR,type='pearson')
Fitted<- fitted(Nest_TR)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

Resid<- resid(poop)
Fitted<- fitted(poop)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

###############################################################################################Relative Food Carrying


######### non transformed data, non homogeneity of variace
FC_TR<-lmer(FoodCarry~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(FC_TR),digits=3)
drop1(FC_TR,test='Chisq')

######### ln transformed data, all good 
FC_TR<-lmer(log(FoodCarry)~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment %>% filter(FoodCarry>0))
print(summary(FC_TR),digits=3)
drop1(FC_TR,test='Chisq')

table(Focal_Treatment %>% filter(FoodCarry>0) %>% select(Treatment))

tidy(FC_TR) %>% 
  write.csv(.,'T_FC_Relative.csv')


FC_TR_PH <- glht(FC_TR, linfct = mcp(Treatment='Tukey'))
summary(FC_TR_PH)


######### log link, not ideal not super homogenous 
Nest_TR<-glmer(Nest~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID), family="gaussian"(link='log'),data=Focal_Treatment %>% filter(Nest>0))
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')


######### Gamma. quite homogenous  
Nest_TR<-glmer(Nest~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=Gamma(link=log),data=Focal_Treatment %>% filter(Nest>0))
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')


######### Rate: terrible
Nest_TR<-lmer(NRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),data=Focal_Treatment)
print(summary(Nest_TR),digits=3)
drop1(Nest_TR,test='Chisq')


######### Rate with beta family
poop<-glmmTMB(NRate~Treatment+ScaledDur+I(ScaledDur^2)+(1|Colony/AnimalID),family=list(family="beta",link="logit"),data=Focal_Treatment %>% filter(Nest>0))



par(mfrow=c(1,1))

Resid<- resid(FC_TR,type='pearson')
Fitted<- fitted(FC_TR)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

Resid<- resid(poop)
Fitted<- fitted(poop)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

#normality of residuals
hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

###############################################################################################Effect on submissive calls and aggression during the focal  
######################################################################################

SI_Pass_All5<-read.csv("SI_Pass_All4.csv") %>% 
  filter(SexPartner=="F") %>% 
filter(Partner_Weight_Min>60) %>% 
  select(Colony,AnimalID,Treatment,DurationSeen,Partner,SexPartner,Partner_Weight_Min,Received,Pass_NB,Bite,Chase,Displace,Overt,Pulltail,Retreat,Bouts_NB,Call_NB,Partner_NB) %>% mutate(SI=Bite+Chase+Displace+Overt+Pulltail,SI_Retreat=Bite+Chase+Displace+Overt+Pulltail+Retreat) %>% 
  mutate(ChaseOvert=Chase+Overt) %>% 
   group_by(Colony,AnimalID,Treatment,DurationSeen,Received) %>%
  mutate(Partner_NB_sum=n()) %>% 
  ungroup() %>% 
  select(-c(Partner,SexPartner,Partner_Weight_Min)) %>%
  group_by(Colony,AnimalID,Treatment,DurationSeen,Partner_NB_sum,Received) %>% 
  summarise_all(sum) %>% 
  ungroup() %>% 
   mutate(Treatment = factor(Treatment, levels = c( "Control", "Tlow","Thigh")))

names(SI_Pass_All5)

Call_Agg_Focal<-SI_Pass_All5 %>% 
  gather(Behaviour,Count,Bite:ChaseOvert)%>% 
  unite(Key,Received,Behaviour,sep="") %>% 
  spread(Key,Count) 
names(Call_Agg_Focal)
View(Call_Agg_Focal)


###############################################################################################overdispersion 
######################################################################################


FUN <- function(fit) {
  #return(fixef(fit))
  x<-resid(fit,type="pearson")
  return(sum(x^2))
}	


od<-function(bootobject){
  biasvals<-bootobject $t0/bootobject[2]$t#$t0 indicates the original value of sum of square of pearson-s residual I had in my minimum model FoodW3 whereas bootobject$t indicates the sum of square of pearson's residual in the simulated data. Thus one ratio or one overdispersion statistic is generated for each simulated dataset
  bias<-mean(biasvals,na.rm=T)
  intervals<-quantile(biasvals,c(0.025,0.975),na.rm=T)
  dat<-c(bias,intervals)
  return(dat)
}

od.point<-function(modelobject){
  x<-sum(resid(modelobject,type="pearson")^2)
  rdf<-summary(modelobject)$AICtab[5]
  return(x/rdf)
}



###############################################################################################Submissive Calls given 
######################################################################################

names(Call_Agg_Focal)


####I have two options, either I go for i) duration and interactions received or ii) interactions receibved and passes. I have decided to go for ii) which is different than in the thesis

MyVar<-c("DurationSeen","Pass_NB","InitiatedSI","ReceivedSI","Partner_NB_sum")
MyVar<-c("DurationSeen","InitiatedSI","ReceivedSI","Partner_NB_sum")
Mypairs(Call_Agg_Focal[,MyVar])

MyVar<-c("DurationSeen","Pass_NB","ReceivedSI")
MyVar<-c("Pass_NB","ReceivedSI")
corvif(Call_Agg_Focal[,MyVar])


####call given, if model simplification social interaction received is not significant any more

library(lme4)

Call_G<-glmer.nb(InitiatedBouts_NB~Treatment+scale(ReceivedSI)+scale(Pass_NB)+(1|Colony/AnimalID),data=Call_Agg_Focal)
summary(Call_G)
drop1(Call_G,test="Chi")

print(summary(Call_G),digits=3)

tidy(Call_G) %>% 
  write.csv(.,'Call_G.csv')

Call_G_PH <- glht(Call_G, linfct = mcp(Treatment='Tukey'))
summary(Call_G_PH)


####Model Validation 
Resid<- resid(Call_G,type='pearson')
Fitted<- fitted(Call_G)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)


Fig4<-ggplot(Call_Agg_Focal,aes(Treatment,InitiatedBouts_NB/DurationSeen))+
  geom_boxplot(aes(fill=Treatment),outlier.shape=NA)+
  geom_point(size=3,alpha=0.3)+
  #geom_jitter(width=0,height=0.2,size=2,col=1,alpha=0.6)+
  geom_line(aes(group=AnimalID),alpha=0.3)+
  scale_x_discrete("",expand=c(0.05,0.05))+
  theme(axis.text.x = element_text(size = 12))+
  scale_y_continuous("Bout of submissive calls given per min",limits=c(0,0.20))+
  theme_classic()+
  scale_fill_manual(values = c("white", "#00CCFF","#FF3300"))+
  guides(fill=FALSE)

####no overdispersion 
od.point(Call_G)

m1boot<-bootMer(Call_G,FUN,100)
library(boot) #required to inspect results
od(m1boot)



###############################################################################################Submissive Calls Received
######################################################################################

####perhaps a zero inflated model would have been more appropriated, I could redo the simulation as in Zuur
table(Call_Agg_Focal$ReceivedBouts_NB)


MyVar<-c("DurationSeen","Pass_NB","Partner_NB_sum","InitiatedBouts_NB")
MyVar<-c("DurationSeen","Partner_NB_sum")
Mypairs(Call_Agg_Focal[,MyVar])
corvif(Call_Agg_Focal[,MyVar])


####call Received
Call_R<-glmer.nb(ReceivedBouts_NB~Treatment+scale(InitiatedSI)+scale(Pass_NB)+(1|Colony/AnimalID),data=Call_Agg_Focal)
summary(Call_R)
drop1(Call_R,test="Chi")


Call_R_PH <- glht(Call_R, linfct = mcp(Treatment='Tukey'))
summary(Call_R_PH)

print(summary(Call_R),digits=3)

tidy(Call_R) %>% 
  write.csv(.,'Call_R.csv')


Resid<- resid(Call_R,type='pearson')
Fitted<- fitted(Call_R)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

####no overdispersion 
od.point(Call_R)

m1boot<-bootMer(Call_R,FUN,100)
library(boot) #required to inspect results
od(m1boot)

###############################################################################################Aggression Given
######################################################################################

MyVar<-c("DurationSeen","Pass_NB","Partner_NB_sum")
MyVar<-c("DurationSeen","Pass_NB")
MyVar<-c("Pass_NB","Partner_NB_sum")
MyVar<-c("DurationSeen","Partner_NB_sum")
Mypairs(Call_Agg_Focal[,MyVar])
corvif(Call_Agg_Focal[,MyVar])


####Aggression given
Agg_G<-glmer.nb(InitiatedSI~Treatment+scale(Pass_NB)+scale(DurationSeen)+(1|Colony/AnimalID),data=Call_Agg_Focal)
summary(Agg_G)
drop1(Agg_G,test="Chi")

print(summary(Agg_G),digits=3)

tidy(Agg_G) %>% 
  write.csv(.,'Agg_G.csv')

Agg_G_PH <- glht(Agg_G, linfct = mcp(Treatment='Tukey'))
summary(Agg_G_PH)


Resid<- resid(Agg_G,type='pearson')
Fitted<- fitted(Agg_G)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

####no overdispersion 
od.point(Agg_G)

m1boot<-bootMer(Agg_G,FUN,100)
library(boot) #required to inspect results
od(m1boot)

Fig5<-ggplot(Call_Agg_Focal,aes(Treatment,InitiatedSI/DurationSeen))+
  geom_boxplot(aes(fill=Treatment),outlier.shape=NA)+
  geom_point(size=3,alpha=0.3)+
  #geom_jitter(width=0,height=0.2,size=2,col=1,alpha=0.6)+
  geom_line(aes(group=AnimalID),alpha=0.3)+
  scale_x_discrete("",expand=c(0.05,0.05))+
  theme(axis.text.x = element_text(size = 12))+
  scale_y_continuous("Aggressive behaviours given per min",limits=c(0,0.15))+
  theme_classic()+
  scale_fill_manual(values = c("white", "#00CCFF","#FF3300"))+
  guides(fill=FALSE)

###############################################################################################Aggression Received
######################################################################################

####Aggression Received
Agg_R<-glmer.nb(ReceivedSI~Treatment+scale(Pass_NB)+scale(DurationSeen)+(1|Colony/AnimalID),data=Call_Agg_Focal)
summary(Agg_R)
drop1(Agg_R,test="Chi")

print(summary(Agg_R),digits=3)

tidy(Agg_R) %>% 
  write.csv(.,'Agg_R.csv')

Agg_R_PH <- glht(Agg_R, linfct = mcp(Treatment='Tukey'))
summary(Agg_R_PH)


Resid<- resid(Agg_R,type='pearson')
Fitted<- fitted(Agg_R)
plot(x = Fitted,
     y = Resid,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 1)

hist(Resid,breaks = 8)
qqnorm(Resid)
qqline(Resid)

####no overdispersion 
od.point(Agg_R)

m1boot<-bootMer(Agg_R,FUN,100)
library(boot) #required to inspect results
od(m1boot)



###############################################################################################Dispersal proxies
######################################################################################


Disp_Proxies<-read.csv('DispersalProxies.csv') %>% 
  mutate(Treatment = factor(Treatment, levels = c( "Control", "Tlow","Thigh"))) 
  
dim(Disp_Proxies)


MyVar<-c("Lat","Unfamiliar_Dur","Activity","Area")
Mypairs(Disp_Proxies[,MyVar])
corvif(Disp_Proxies[,MyVar])


########Latency
Latency<-lmer(log(Lat)~Treatment+Treatment_NB+(1|Colony/AnimalID),data=Disp_Proxies)
summary(Latency)
drop1(Latency,test="Chisq")
print(summary(Latency),digits=3)

tidy(Latency) %>% 
  write.csv(.,'T_Latency.csv')

Latency_PH <- glht(Latency, linfct = mcp(Treatment='Tukey'))
summary(Latency_PH)

Valid(Latency)

########duration in unfamiliar location
Unfamiliar<-lmer(log(Unfamiliar_Dur)~Treatment+Treatment_NB+(1|Colony/AnimalID),data=Disp_Proxies %>% filter (Unfamiliar_Dur>0))
summary(Unfamiliar)
drop1(Unfamiliar,test="Chisq")
print(summary(Unfamiliar),digits=3)

table(Disp_Proxies %>% filter (Unfamiliar_Dur>0) %>% select(Treatment))

tidy(Unfamiliar) %>% 
  write.csv(.,'T_Unfamiliar.csv')

Unfamiliar_PH <- glht(Unfamiliar, linfct = mcp(Treatment='Tukey'))
summary(Unfamiliar_PH)

Valid(Unfamiliar)

Fig6<-ggplot(Disp_Proxies,aes(Treatment,Unfamiliar_Dur))+
  geom_boxplot(aes(fill=Treatment),outlier.shape=NA)+
  geom_point(size=3,alpha=0.3)+
  #geom_jitter(width=0,height=0.2,size=2,col=1,alpha=0.6)+
  geom_line(aes(group=AnimalID),alpha=0.3)+
  scale_x_discrete("",expand=c(0.05,0.05))+
  theme(axis.text.x = element_text(size = 12))+
  scale_y_continuous("Total time spent in unfamiliar areas in sec",limits=c(0,600))+
  theme_classic()+
  scale_fill_manual(values = c("white", "#00CCFF","#FF3300"))+
  guides(fill=FALSE)


########Total nb of area explored 
Activity<-glmer.nb(Activity~Treatment+Treatment_NB+(1|Colony/AnimalID),data=Disp_Proxies)
summary(Activity)
drop1(Activity,test="Chisq")
print(summary(Activity),digits=3)

tidy(Activity) %>% 
  write.csv(.,'T_Activity.csv')

Activity_PH <- glht(Activity, linfct = mcp(Treatment='Tukey'))
summary(Activity_PH)

Valid(Activity)

####no overdispersion 
od.point(Activity)

m1boot<-bootMer(Activity,FUN,100)
library(boot) #required to inspect results
od(m1boot)

##############################################################################################Breeding female encounter
######################################################################################


##############################################################################################Non familiar breeding female encounter
######################################################################################




