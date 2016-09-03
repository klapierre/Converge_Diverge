library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

#kim's
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

###read in data

#experiment information
expInfo <- read.csv('ExperimentInformation_Mar2016.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep='::'))%>%
  select(-X)%>%
  filter(treatment_year!=0)

#diversity data
div <- merge(read.csv('DiversityMetrics_Sept2016.csv'), expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  select(-X)%>%
  filter(treatment_year!=0)

anpp<-read.csv("ANPP_March2016.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)

###calculate change in dispersion, H, S, and evenness

#subset out controls and treatments
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, dispersion, H, S, SimpEven)
  names(divControls)[names(divControls)=='dispersion'] <- 'ctl_dispersion'
  names(divControls)[names(divControls)=='H'] <- 'ctl_H'
  names(divControls)[names(divControls)=='S'] <- 'ctl_S'
  names(divControls)[names(divControls)=='SimpEven'] <- 'ctl_SimpEven'
divTrt <- div%>%
  #removing treatments that were low levels of resource manipulation, not manipulating any resource, pulses
  #plus filtering to get only treatments
  filter(max_trt==1, resource_mani==1, pulse==0, plot_mani>0)

#merge controls and treatments
divCompare <- merge(divControls, divTrt, by=c('exp_year'))%>%
#calculate change in disperion, H, S, and evenness
  mutate(dispersion_change=dispersion-ctl_dispersion, 
         H_change=H-ctl_H, 
         S_PC=(S-ctl_S)/ctl_S, 
         SimpEven_change=SimpEven-ctl_SimpEven)%>%
  select(exp_year, treatment_year, treatment, plot_mani, mean_change, dispersion_change, H_change,  SimpEven_change, S_PC, site_code, project_name, community_type, calendar_year)
# 
##comparing change vs percent change
# d1<-qplot(dispersion_PC, data=divCompare, geom="histogram")+
#   ggtitle("dispersion percent change")

theme_set(theme_bw(16))
d2<-qplot(dispersion_change, data=divCompare, geom="histogram")+
  ggtitle("Within Treatment Change")+
  xlab("Trt Disp - Cont Disp")+
  geom_vline(xintercept = 0, size=2)

m<-qplot(mean_change, data=divCompare, geom="histogram")+
  ggtitle("Among Treatment Change")+
  xlab(" Distance between Centriods")+
  geom_vline(xintercept = 0, size=2)

s1<-qplot(S_PC, data=divCompare, geom="histogram")+
  ggtitle("Richness Percent Change")+
  xlab("Percent Change in Richness")+
  geom_vline(xintercept = 0, size=2)
# s2<-qplot(S_change, data=divCompare, geom="histogram")+
#   ggtitle("richness change")

# e1<-qplot(SimpEven_PC, data=divCompare, geom="histogram")+
#   ggtitle("even percent change")
e2<-qplot(SimpEven_change, data=divCompare, geom="histogram")+
  ggtitle("Evenness Change")+
  xlab("Trt Evenness - Cont Evenness")+
  geom_vline(xintercept = 0, size=2)

grid.arrange( m, d2,s1, e2, ncol=2)

# #without cdr
# divCompareTest <- divCompare%>%
#   filter(site_code!='CDR')
# 
# theme_set(theme_bw(16))
# d2<-qplot(dispersion_change, data=divCompareTest, geom="histogram")+
#   ggtitle("Within Treatment Change")+
#   xlab("Trt Disp - Cont Disp")+
#   geom_vline(xintercept = 0, size=2)
# 
# m<-qplot(mean_change, data=divCompareTest, geom="histogram")+
#   ggtitle("Among Treatment Change")+
#   xlab(" Distance between Centriods")+
#   geom_vline(xintercept = 0, size=2)
# 
# s1<-qplot(S_PC, data=divCompareTest, geom="histogram")+
#   ggtitle("Richness Percent Change")+
#   xlab("Percent Change in Richness")+
#   geom_vline(xintercept = 0, size=2)
# # s2<-qplot(S_change, data=divCompare, geom="histogram")+
# #   ggtitle("richness change")
# 
# # e1<-qplot(SimpEven_PC, data=divCompare, geom="histogram")+
# #   ggtitle("even percent change")
# e2<-qplot(SimpEven_change, data=divCompareTest, geom="histogram")+
#   ggtitle("Evenness Change")+
#   xlab("Trt Evenness - Cont Evenness")+
#   geom_vline(xintercept = 0, size=2)
# 
# grid.arrange( m, d2,s1, e2, ncol=2)


###merging with experiment (treatment) information
SiteExp<-read.csv("SiteExperimentDetails_March2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))

#full dataset
write.csv(ForAnalysis, "ForBayesianAnalysis_Sept2016.csv")

#9 yr or less
ForAnalysis9yr <- ForAnalysis%>%
  mutate(treatment_year=ifelse(project_name=='TMECE'&treatment_year==11, 1, ifelse(project_name=='TMECE'&treatment_year==12, 2, ifelse(project_name=='TMECE'&treatment_year==13, 3, ifelse(project_name=='TMECE'&treatment_year==14, 4, ifelse(project_name=='TMECE'&treatment_year==15, 5, ifelse(project_name=='TMECE'&treatment_year==16, 6, ifelse(project_name=='TMECE'&treatment_year==17, 7, ifelse(project_name=='TMECE'&treatment_year==18, 8, ifelse(project_name=='TMECE'&treatment_year==19, 9, ifelse(project_name=='TMECE'&treatment_year==20, 10, ifelse(project_name=='TMECE'&treatment_year==21, 11, ifelse(project_name=='TMECE'&treatment_year==22, 12, ifelse(project_name=='TMECE'&treatment_year==23, 13, ifelse(project_name=='TMECE'&treatment_year==24, 14, ifelse(project_name=='TMECE'&treatment_year==25, 15, ifelse(project_name=='TMECE'&treatment_year==26, 16, ifelse(project_name=='TMECE'&treatment_year==27, 17, treatment_year))))))))))))))))))%>%
  filter(treatment_year<10)
write.csv(ForAnalysis9yr, "ForBayesianAnalysis_9yr_Sept2016.csv")

#10+ year datasets (all years)
ForAnalysis10yr <- ForAnalysis%>%
  filter(experiment_length>9)
write.csv(ForAnalysis10yr, "ForBayesianAnalysis_10yr_Sept2016.csv")

#absolute value
ForAnalysisAbsValue <- ForAnalysis%>%
  mutate(mean_change=abs(mean_change), dispersion_change=abs(dispersion_change), H_change=abs(H_change), SimpEven_change=abs(SimpEven_change), S_PC=abs(S_PC))
write.csv(ForAnalysisAbsValue, "ForBayesianAnalysis_abs value_9yr_Sept2016.csv")








##doing the same thing for anpp
anppMeans<-anpp%>%
  tbl_df()%>%
  group_by(site_code, project_name, calendar_year, treatment, community_type)%>%
  summarize(anpp=mean(anpp))

anppMeans2<-merge(anppMeans, expInfo, by=c("site_code","project_name","community_type","treatment", "calendar_year"))#rows are dropped b/c species were not recorded those years

anppControls <- subset(anppMeans2, subset=(plot_mani==0))%>%
  select(exp_year, anpp)
names(anppControls)[names(anppControls)=='anpp'] <- 'ctl_anpp'
anppTrt <- subset(anppMeans2, subset=(plot_mani!=0))

#merge controls and treatments
anppCompare <- merge(anppControls, anppTrt, by=c('exp_year'))%>%
  mutate(anpp_PC=(anpp-ctl_anpp)/ctl_anpp)%>%
  select(exp_year, treatment, plot_mani, anpp_PC)

# a1<-qplot(anpp_PC, data=anppCompare, geom="histogram")+
#   ggtitle("anpp percent change")
# a2<-qplot(anpp_change, data=anppCompare, geom="histogram")+
#   ggtitle("anpp change")
# 
# grid.arrange(d1, d2, s1, s2, e1, e2, a1, a2, ncol=2)

###merging with experiment (treatment) information
anppCompareExp <- merge(anppCompare, expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  #removing treatments that were pulses, did not directly manipulate a resource, or had ceased and pre-treatment data
  filter(pulse==0, resource_mani==1, max_trt==1, treatment_year>0)%>%
  select(exp_year, treatment, plot_mani, anpp_PC, site_code, project_name, community_type, calendar_year, treatment_year)

ForANPPAnalysis<-merge(anppCompareExp, SiteExp, by=c("site_code","project_name","community_type"))

qplot(anpp_PC, data=anppCompareExp, geom="histogram")+
  xlab("ANPP Percent Change")+
  geom_vline(xintercept = 0, size=2)

test<-ForANPPAnalysis%>%
  select(site_code, project_name, community_type)%>%
  unique()

write.csv(ForANPPAnalysis, "ForBayesianAnalysisANPP_March2016b.csv")
