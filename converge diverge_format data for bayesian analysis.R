library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

###read in data

#experiment information
expInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X)

anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)

#diversity data
div <- read.csv('DiversityMetrics_May2017.csv')%>%
  separate(exp_year, c('site_code', 'project_name','community_type', 'calendar_year'), sep='::')%>%
  mutate(calendar_year=as.integer(calendar_year))%>%
  left_join(expInfo)%>%
  select(-X)%>%
  filter(treatment_year!=0)

#anpp data
anpp<-read.csv("ANPP_Oct2017.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)

#appearance/disappearance data
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

# turnover <- read.csv('appear_disappear_Mar2017.csv')%>%
#   separate(col=exp_code, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::')%>%
#   left_join(expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment', 'calendar_year'), all=T)%>%
#   filter(pulse==0&project_name!='e001'&project_name!='e002')
# turnoverCDRe001 <- turnover%>%filter(project_name=='e001'&treatment!=2&treatment!=3&treatment!=4&treatment!=5&treatment!=7)
# turnoverCDRe002 <- turnover%>%filter(project_name=='e002'&treatment!='2_f_u_n'&treatment=='3_f_u_n'&treatment=='4_f_u_n'&treatment=='5_f_u_n'&treatment=='7_f_u_n')
# turnoverAll<-rbind(turnover, turnoverCDRe001, turnoverCDRe002)%>%
#   left_join(SiteExp, by=c('site_code', 'project_name', 'community_type'))
# 
# #full dataset
# # write.csv(turnoverAll, "turnover_all years_Mar2017.csv")
# 
# #9 yr or less
# turnover9yr <- turnoverAll%>%
#   filter(treatment_year<10)
# 
# ##18% of our our data is from CDR and KNZ
# 
# # write.csv(turnover9yr, "turnover_9yr_May2017.csv")
# 
# ggplot(turnover9yr, aes(x=treatment_year, y=appearance/total, color=plot_mani)) +
#   geom_point()
# ggplot(turnover9yr, aes(x=treatment_year, y=disappearance/total, color=plot_mani)) +
#   geom_point()



###calculate change in dispersion, H, S, and evenness

#subset out controls and treatments
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, dispersion, H, S, SimpEven, calendar_year, treatment_year)
  names(divControls)[names(divControls)=='dispersion'] <- 'ctl_dispersion'
  names(divControls)[names(divControls)=='H'] <- 'ctl_H'
  names(divControls)[names(divControls)=='S'] <- 'ctl_S'
  names(divControls)[names(divControls)=='SimpEven'] <- 'ctl_SimpEven'
divTrt1 <- div%>%
  #removing treatments that were low levels of resource manipulation, not manipulating any resource, pulses
  #plus filtering to get only treatments
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")

#removing a subset of CDR treatments to prevent the majority of data being from CDR; keeping lowest, highest, and 10 gm-2 (level most comparable to other studies)
divCDRe001<-div%>%
  filter(site_code=="CDR"&treatment==1|treatment==6|treatment==8|treatment==9,plot_mani>0)
divCDRe002<-div%>%
  filter(site_code=="CDR"&treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n',plot_mani>0)
  

divTrt<-rbind(divTrt1, divCDRe002, divCDRe001)

##16% of our data is from CDR and 10% is from KNZ
#merge controls and treatments
divCompare <- divControls%>%
  left_join(divTrt)%>%
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
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))

#full dataset
# write.csv(ForAnalysis, "ForBayesianAnalysis_May2017.csv")

#8 yr or less
ForAnalysis8yr <- ForAnalysis%>%
  filter(treatment_year<9)

##18% of our our data is from CDR and KNZ

# write.csv(ForAnalysis8yr, "ForBayesianAnalysis_8yr_May2017.csv")


#Plot of 8 year data used in paper.
theme_set(theme_bw(16))
d2<-qplot(dispersion_change, data=ForAnalysis8yr, geom="histogram")+
  ggtitle("Within Treatment Change")+
  xlab("Trt Disp - Cont Disp")+
  geom_vline(xintercept = 0, size=2)

m<-qplot(mean_change, data=ForAnalysis8yr, geom="histogram")+
  ggtitle("Among Treatment Change")+
  xlab(" Distance between Centriods")+
  geom_vline(xintercept = 0, size=2)

s1<-qplot(S_PC, data=ForAnalysis8yr, geom="histogram")+
  ggtitle("Richness Percent Change")+
  xlab("Percent Change in Richness")+
  geom_vline(xintercept = 0, size=2)
# s2<-qplot(S_change, data=divCompare, geom="histogram")+
#   ggtitle("richness change")

# e1<-qplot(SimpEven_PC, data=divCompare, geom="histogram")+
#   ggtitle("even percent change")
e2<-qplot(SimpEven_change, data=ForAnalysis8yr, geom="histogram")+
  ggtitle("Evenness Change")+
  xlab("Trt Evenness - Cont Evenness")+
  geom_vline(xintercept = 0, size=2)

grid.arrange( m, d2,s1, e2, ncol=2)


# #look at pre-trt data
# pre <- ForAnalysis8yr%>%
#   filter(site_code=='KUFS'|site_code=='GVN'|site_code=='dcgs'|site_code=='JSP'|site_code=='BAY'|site_code=='CEH'|site_code=='PIE'|site_code=='KAEFS'|project_name=='pplots'|project_name=='snow')
# 
# theme_set(theme_bw(16))
# d2<-qplot(dispersion_change, data=subset(pre, treatment_year==1), geom="histogram")+
#   ggtitle("Within Treatment Change")+
#   xlab("Trt Disp - Cont Disp")+
#   geom_vline(xintercept = 0, size=2)
# 
# m<-qplot(mean_change, data=subset(pre, treatment_year==1), geom="histogram")+
#   ggtitle("Among Treatment Change")+
#   xlab(" Distance between Centriods")+
#   geom_vline(xintercept = 0, size=2)
# 
# s1<-qplot(S_PC, data=subset(pre, treatment_year==1), geom="histogram")+
#   ggtitle("Richness Percent Change")+
#   xlab("Percent Change in Richness")+
#   geom_vline(xintercept = 0, size=2)
# # s2<-qplot(S_change, data=divCompare, geom="histogram")+
# #   ggtitle("richness change")
# 
# # e1<-qplot(SimpEven_PC, data=divCompare, geom="histogram")+
# #   ggtitle("even percent change")
# e2<-qplot(SimpEven_change, data=subset(pre, treatment_year==1), geom="histogram")+
#   ggtitle("Evenness Change")+
#   xlab("Trt Evenness - Cont Evenness")+
#   geom_vline(xintercept = 0, size=2)
# 
# grid.arrange( m, d2,s1, e2, ncol=2)
# 
# ggplot(data=pre, aes(x=treatment_year, y=mean_change)) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=site_code))

#9+ year datasets (all years)
ForAnalysis9yr <- ForAnalysis%>%
  filter(experiment_length>8)
# write.csv(ForAnalysis9yr, "ForBayesianAnalysis_9plusyr_May2017.csv")

#absolute value
ForAnalysisAbsValue <- ForAnalysis9yr%>%
  mutate(mean_change=abs(mean_change), dispersion_change=abs(dispersion_change), H_change=abs(H_change), SimpEven_change=abs(SimpEven_change), S_PC=abs(S_PC))
# write.csv(ForAnalysisAbsValue, "ForBayesianAnalysis_abs value_9yr_Dec2016.csv")



###getting the treatment interaction types
trtType <- ForAnalysis%>%
  left_join(expInfo)%>%
  #create drought and irrigation categories
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(plot_mani==2&n>0&drought<0, 'N+drought', ifelse(plot_mani==2&n>0&irrigation>0, 'N+irr', ifelse(plot_mani==2&n>0&p>0, 'N+P', ifelse(plot_mani==2&p>0&k>0, 'P+K', ifelse(plot_mani==2&n>0&CO2>0, 'N+CO2', ifelse(plot_mani==2&CO2>0&irrigation>0, 'CO2+irr', ifelse(plot_mani==3&n>0&p>0&k>0, 'N+P+K', ifelse(plot_mani==3&n>0&CO2>0&irrigation>0, 'N+CO2+irr', ifelse(plot_mani==4&n>0&p>0&k>0&irrigation>0, 'N+P+K+irr', ifelse(plot_mani==1&n>0, 'N', ifelse(plot_mani==1&p>0, 'P', ifelse(plot_mani==1&irrigation>0, 'irr', ifelse(plot_mani==1&drought<0, 'drought', ifelse(plot_mani==1&CO2>0, 'CO2', ifelse(n==0&drought==0&irrigation==0&p==0&k==0&CO2==0, 'other', 'resource+other'))))))))))))))))%>%
  #drop experiments that we can't run the analyses on
  filter(plot_mani<6, treatment_year!=0, anpp!='NA')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, dispersion_change, SimpEven_change, S_PC, experiment_length, resource_other)%>%
  #keep final year only
  group_by(site_code, project_name, community_type, treatment)%>%
  filter(treatment_year==max(treatment_year))%>%
  ungroup()

# write.csv(trtType, 'treatment interactions_11152017.csv')

###getting resource*non-resource interactions
trtTypeRes <- ForAnalysis%>%
  left_join(expInfo)%>%
  #create drought and irrigation categories
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(plot_mani==2&n>0&drought<0, 'N+drought', ifelse(plot_mani==2&n>0&irrigation>0, 'N+irr', ifelse(plot_mani==2&n>0&p>0, 'N+P', ifelse(plot_mani==2&p>0&k>0, 'P+K', ifelse(plot_mani==2&n>0&CO2>0, 'N+CO2', ifelse(plot_mani==2&CO2>0&irrigation>0, 'CO2+irr', ifelse(plot_mani==3&n>0&p>0&k>0, 'N+P+K', ifelse(plot_mani==3&n>0&CO2>0&irrigation>0, 'N+CO2+irr', ifelse(plot_mani==4&n>0&p>0&k>0&irrigation>0, 'N+P+K+irr', ifelse(plot_mani==1&n>0, 'N', ifelse(plot_mani==1&p>0, 'P', ifelse(plot_mani==1&irrigation>0, 'irr', ifelse(plot_mani==1&drought<0, 'drought', ifelse(plot_mani==1&CO2>0, 'CO2', ifelse(n==0&drought==0&irrigation==0&p==0&k==0&CO2==0, 'other', 'resource+other'))))))))))))))))%>%
  #create resource*other trt column
  mutate(resource_other=ifelse(resource_mani!=0&other_trt!=0, 'R*other', ifelse(resource_mani!=0&mow_clip==1, 'R*mow_clip', ifelse(resource_mani==1&burn==1, 'R*burn', ifelse(resource_mani!=0&herb_removal==1, 'R*herbrem', ifelse(resource_mani!=0&temp>0, 'R*temp', ifelse(resource_mani!=0&plant_mani==1, 'R*plant_mani', ifelse(resource_mani==0, 'non-resource', ifelse(resource_mani!=0&plot_mani>1, 'R*R', ifelse(project_name=='e001', 'CDR', ifelse(project_name=='e002', 'CDR', 'single-resource')))))))))))%>%
  #create column of non-resource manipulations; NOTE: this does not account for trts with multiple at the same time, they are overwritten by last one, so only use this column when dropping plot mani>1
  mutate(nonresource=ifelse(mow_clip==1, 'mow_clip', ifelse(burn==1, 'burn', ifelse(herb_removal==1, 'herbrem', ifelse(temp>0, 'temp', ifelse(plant_mani==1, 'plant_mani', 'other'))))))%>%
  #drop experiments that we can't run the analyses on
  filter(plot_mani<6, treatment_year!=0, anpp!='NA')%>%
  #drop precip variability treatments
  filter(other_trt!='reduced precip variability'&other_trt!='increased precip variability'&other_trt!='increase winter precip, decrease summer precip'&other_trt!='decrease winter precip, increase summer precip')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, resource_other, nonresource, plot_mani, resource_mani, mean_change, dispersion_change, SimpEven_change, S_PC, experiment_length)

#only run to the select line in making trtTypeRes to allow next lines to work
temp <- trtTypeRes%>%
  select(site_code, project_name, treatment, trt_type, resource_other, nonresource, other_trt, mow_clip, burn, herb_removal, temp, plant_mani, plot_mani)%>%
  group_by(site_code, project_name, treatment, trt_type, resource_other, nonresource)%>%
  unique()%>%
  ungroup()

# write.csv(temp, 'treatment_resource_nonresource.csv')

#mean change
#notes: remove KBS because tilling has a big effect and it is the only tilled experiment
singleResourceFig <- ggplot(data=subset(trtTypeRes, plot_mani==1&resource_mani==1&site_code!='KBS'&treatment_year<9), aes(x=treatment_year, y=mean_change)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('Mean Change') +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

singleNonresourceFig <- ggplot(data=subset(trtTypeRes, plot_mani==1&resource_mani==0&site_code!='KBS'&treatment_year<9), aes(x=treatment_year, y=mean_change)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=nonresource)) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('Mean Change') +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

resourceNonresourceFig <- ggplot(data=subset(trtTypeRes, plot_mani>1&resource_mani!=0&site_code!='KBS'&treatment_year<9), aes(x=treatment_year, y=mean_change)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=resource_other)) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('Mean Change') +
  scale_y_continuous(limits=c(0,1)) +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

pushViewport(viewport(layout=grid.layout(1,3)))
print(singleResourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(singleNonresourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(resourceNonresourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))

#richness
#notes: remove KBS because tilling has a big effect and it is the only tilled experiment
singleResourceFig <- ggplot(data=subset(trtTypeRes, plot_mani==1&resource_mani==1&site_code!='KBS'&treatment_year<9), aes(x=treatment_year, y=S_PC)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('Richness Change') +
  scale_y_continuous(limits=c(-1,1.3)) +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

singleNonresourceFig <- ggplot(data=subset(trtTypeRes, plot_mani==1&resource_mani==0&site_code!='KBS'&treatment_year<9), aes(x=treatment_year, y=S_PC)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=nonresource)) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('Richness Change') +
  scale_y_continuous(limits=c(-1,1.3)) +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

resourceNonresourceFig <- ggplot(data=subset(trtTypeRes, plot_mani>1&resource_mani!=0&site_code!='KBS'&treatment_year<9), aes(x=treatment_year, y=S_PC)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=resource_other)) +
  # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
  xlab('Treatment Year') + ylab('Richness Change') +
  scale_y_continuous(limits=c(-1,1.3)) +
  theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))

pushViewport(viewport(layout=grid.layout(1,3)))
print(singleResourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(singleNonresourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(resourceNonresourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))


# #dispersion
# #notes: remove KBS because tilling has a big effect and it is the only tilled experiment
# singleResourceFig <- ggplot(data=subset(trtTypeRes, plot_mani==1&resource_mani==1&site_code!='KBS'), aes(x=treatment_year, y=dispersion_change)) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=trt_type)) +
#   # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
#   xlab('Treatment Year') + ylab('Dispersion Change') +
#   scale_y_continuous(limits=c(-0.5,0.4)) +
#   theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))
# 
# singleNonresourceFig <- ggplot(data=subset(trtTypeRes, plot_mani==1&resource_mani==0&site_code!='KBS'), aes(x=treatment_year, y=dispersion_change)) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=nonresource)) +
#   # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
#   xlab('Treatment Year') + ylab('Dispersion Change') +
#   scale_y_continuous(limits=c(-0.5,0.4)) +
#   theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))
# 
# resourceNonresourceFig <- ggplot(data=subset(trtTypeRes, plot_mani>1&resource_mani!=0&site_code!='KBS'), aes(x=treatment_year, y=dispersion_change)) +
#   geom_smooth(method='lm', formula=y~x+I(x^2), se=F, aes(color=resource_other)) +
#   # geom_smooth(method='lm', formula=y~x+I(x^2), size=3, color='black') +
#   xlab('Treatment Year') + ylab('Dispersion Change') +
#   scale_y_continuous(limits=c(-0.5,0.4)) +
#   theme(legend.position=c(0.05,0.95), legend.justification=c(0,1))
# 
# pushViewport(viewport(layout=grid.layout(1,3)))
# print(singleResourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(singleNonresourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(resourceNonresourceFig, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))

# ANPP data ---------------------------------------------------------------

##doing the same thing for anpp
anppMeans<-anpp%>%
  tbl_df()%>%
  group_by(site_code, project_name, calendar_year, treatment, community_type)%>%
  summarize(anpp=mean(anpp))

anppMeans2<-merge(anppMeans, expInfo, by=c("site_code","project_name","community_type","treatment"))#rows are dropped b/c species were not recorded those years

anppControls <- subset(anppMeans2, subset=(plot_mani==0))%>%
  select(exp_year, calendar_year, anpp)
names(anppControls)[names(anppControls)=='anpp'] <- 'ctl_anpp'
anppTrt <- subset(anppMeans2, subset=(plot_mani!=0))

#merge controls and treatments
anppCompare <- merge(anppControls, anppTrt, by=c('exp_year', 'calendar_year'))%>%
  mutate(anpp_PC=(anpp-ctl_anpp)/ctl_anpp)%>%
  select(exp_year, calendar_year, treatment, plot_mani, anpp_PC)

# a1<-qplot(anpp_PC, data=anppCompare, geom="histogram")+
#   ggtitle("anpp percent change")
# a2<-qplot(anpp_change, data=anppCompare, geom="histogram")+
#   ggtitle("anpp change")
# 
# grid.arrange(d1, d2, s1, s2, e1, e2, a1, a2, ncol=2)

###merging with experiment (treatment) information
anppCompareExp1 <- merge(anppCompare, expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  #removing treatments that were pulses, did not directly manipulate a resource, or had ceased and pre-treatment data
  filter(pulse==0, site_code!="CDR")%>%
  select(exp_year, treatment, plot_mani, anpp_PC, site_code, project_name, community_type, calendar_year)

anppcdre001<-merge(anppCompare, expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  filter(site_code=="CDR"&treatment==1|treatment==6|treatment==8|treatment==9,plot_mani>0)%>%
  select(exp_year, treatment, plot_mani, anpp_PC, site_code, project_name, community_type, calendar_year)
anppcdre002<-merge(anppCompare, expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  filter(site_code=="CDR"&treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n',plot_mani>0)%>%
  select(exp_year, treatment, plot_mani, anpp_PC, site_code, project_name, community_type, calendar_year)

anppCompareExp<-rbind(anppCompareExp1, anppcdre002, anppcdre001)

ForANPPAnalysis<-merge(anppCompareExp, SiteExp, by=c("site_code","project_name","community_type"))

write.csv(ForANPPAnalysis, "ForBayesianAnalysisANPP_Oct2017.csv")
  
ForANPPAnalysis9yr<-ForANPPAnalysis%>%
  filter(treatment_year<10)
write.csv(ForANPPAnalysis9yr, "ForBayesianAnalysisANPP_9yr_Dec2016.csv")


###looking at stability

#calculate CV for each treatment year combo for spatial and temporal. do not do varience, too variable.

anpp_spatial<-anpp%>%
  group_by(site_code, project_name, community_type, treatment_year, calendar_year, treatment)%>%
  summarize(anpp_sp_sd=sd(anpp, na.rm=T),
            anpp_sp_mean=mean(anpp, na.rm=T),
            anpp_sp_cv=(anpp_sp_sd/anpp_sp_mean)*100)%>%
  select(-anpp_sp_sd, -anpp_sp_mean)%>%
  filter(treatment_year!=0)

anpp_spatial_details<-merge(anpp_spatial, anpp_expInfo, by=c("site_code", "project_name", "community_type","treatment"))

anpp_spatial_forAnalysis<-merge(anpp_spatial_details, SiteExp, by=c("site_code", "project_name", "community_type"))%>%
  filter(pulse==0)

write.csv(anpp_spatial_forAnalysis, "~/Dropbox/converge_diverge/datasets/LongForm/ANPP_Spatail_ForAnalysis.csv")

#checking how long the dataset it. drop KNZ_GFP and NANT_wet
todrop<-anpp_temp_cv<-anpp%>%
  select(site_code, project_name, community_type, calendar_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(n=length(calendar_year))%>%
  mutate(drop=ifelse(n<3,1,0))

anpp_sub<-merge(todrop, anpp, by=c("site_code","project_name","community_type"))%>%
  filter(drop!=1)

anpp_temp_cv<-anpp_sub%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_cv=(anpp_temp_sd/anpp_temp_mean)*100)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(anpp_temp_cv=mean(anpp_temp_cv, na.rm=T))

anpp_temp_details<-merge(anpp_temp_cv, anpp_expInfo, by=c("site_code", "project_name", "community_type","treatment"))

anpp_temp_forAnalysis<-merge(anpp_temp_details, SiteExp, by=c("site_code", "project_name", "community_type"))%>%
  filter(pulse==0)

write.csv(anpp_temp_forAnalysis, "~/Dropbox/converge_diverge/datasets/LongForm/ANPP_Temporal_ForAnalysis.csv")


# ForANPPAnalysis9yr<-ForANPPAnalysis%>%
#   filter(treatment_year<10)
# write.csv(ForANPPAnalysis9yr, "ForBayesianAnalysisANPP_9yr_Dec2016.csv")
# 
# qplot(anpp_PC, data=anppCompareExp, geom="histogram")+
#   xlab("ANPP Percent Change")+
#   geom_vline(xintercept = 0, size=2)
# 
# test<-ForANPPAnalysis%>%
#   select(site_code, project_name, community_type)%>%
#   unique()
# 
# 

###ANPP to MAP
test<-ForAnalysis%>%
  select(site_code, MAP, anpp)%>%
  unique()

plot(MAP, anpp, data=test)
