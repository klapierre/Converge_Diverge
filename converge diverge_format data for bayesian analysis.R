library(tidyr)
library(dplyr)

#kim's
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

###read in data

#experiment information
expInfo <- read.csv('ExperimentInformation_Feb2016.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep='::'))

#diversity data
div <- read.csv('DiversityMetrics_Feb2016.csv')

anpp<-read.csv("ANPP_Feb2016.csv")%>%
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
divTrt <- subset(div, subset=(plot_mani!=0))

#merge controls and treatments
divCompare <- merge(divControls, divTrt, by=c('exp_year'))%>%
#calculate change in disperion, H, S, and evenness
  mutate(dispersion_change=dispersion-ctl_dispersion, H_change=H-ctl_H, S_change=S-ctl_S, SimpEven_change=SimpEven-ctl_SimpEven)%>%
  select(exp_year, treatment, plot_mani, mean_change, dispersion_change, H_change, S_change, SimpEven_change)

###merging with experiment (treatment) information
divCompareExp <- merge(divCompare, expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  #removing treatments that were pulses, did not directly manipulate a resource, or had ceased and pre-treatment data
  filter(pulse==0, resource_mani==1, treatment_year>0)%>%
  select(exp_year, treatment, plot_mani, mean_change, dispersion_change, H_change, S_change, SimpEven_change, site_code, project_name, community_type, calendar_year)

SiteExp<-read.csv("SiteExperimentDetails_Feb2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompareExp, SiteExp, by=c("site_code","project_name","community_type"))

write.csv(ForAnalysis, "ForBayesianAnalysis_Feb2016.csv")


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
  #calculate change in disperion, H, S, and evenness
  mutate(anpp_change=anpp-ctl_anpp)%>%
  select(exp_year, treatment, plot_mani, anpp_change)

###merging with experiment (treatment) information
anppCompareExp <- merge(anppCompare, expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  #removing treatments that were pulses, did not directly manipulate a resource, or had ceased and pre-treatment data
  filter(pulse==0, resource_mani==1, treatment_year>0)%>%
  select(exp_year, treatment, plot_mani, anpp_change, site_code, project_name, community_type, calendar_year)

ForANPPAnalysis<-merge(anppCompareExp, SiteExp, by=c("site_code","project_name","community_type"))

test<-ForANPPAnalysis%>%
  select(site_code, project_name, community_type)%>%
  unique()

write.csv(ForANPPAnalysis, "ForBayesianAnalysisANPP_Feb2016.csv")
