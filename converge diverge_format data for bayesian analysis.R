library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

###read in data

#experiment information
expInfo <- read.csv('ExperimentInformation_ANPP_Oct2017.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X)

anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)

#diversity data
div <- merge(read.csv('DiversityMetrics_May2017.csv'), expInfo, by=c('exp_year', 'treatment', 'plot_mani'))%>%
  select(-X)%>%
  filter(treatment_year!=0)

#anpp data
<<<<<<< HEAD
anpp<-read.csv("ANPP_Oct2017.csv")%>%
=======
anpp<-read.csv("ANPP_OCT2017.csv")%>%
>>>>>>> 455efdf31caabc67e5a48e71fa15195d60f20a05
  select(-X)%>%
  filter(treatment_year!=0)

#appearance/disappearance data
SiteExp<-read.csv("SiteExperimentDetails_March2016.csv")%>%
  select(-X)

turnover <- read.csv('appear_disappear_Mar2017.csv')%>%
  separate(col=exp_code, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::')%>%
  left_join(expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment', 'calendar_year'), all=T)%>%
  filter(pulse==0&project_name!='e001'&project_name!='e002')
turnoverCDRe001 <- turnover%>%filter(project_name=='e001'&treatment!=2&treatment!=3&treatment!=4&treatment!=5&treatment!=7)
turnoverCDRe002 <- turnover%>%filter(project_name=='e002'&treatment!='2_f_u_n'&treatment=='3_f_u_n'&treatment=='4_f_u_n'&treatment=='5_f_u_n'&treatment=='7_f_u_n')
turnoverAll<-rbind(turnover, turnoverCDRe001, turnoverCDRe002)%>%
  left_join(SiteExp, by=c('site_code', 'project_name', 'community_type'))

#full dataset
write.csv(turnoverAll, "turnover_all years_Mar2017.csv")

#9 yr or less
turnover9yr <- turnoverAll%>%
  filter(treatment_year<10)

##18% of our our data is from CDR and KNZ

write.csv(turnover9yr, "turnover_9yr_May2017.csv")

ggplot(turnover9yr, aes(x=treatment_year, y=appearance/total, color=plot_mani)) +
  geom_point()
ggplot(turnover9yr, aes(x=treatment_year, y=disappearance/total, color=plot_mani)) +
  geom_point()



###calculate change in dispersion, H, S, and evenness

#subset out controls and treatments
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, dispersion, H, S, SimpEven)
  names(divControls)[names(divControls)=='dispersion'] <- 'ctl_dispersion'
  names(divControls)[names(divControls)=='H'] <- 'ctl_H'
  names(divControls)[names(divControls)=='S'] <- 'ctl_S'
  names(divControls)[names(divControls)=='SimpEven'] <- 'ctl_SimpEven'
divTrt1 <- div%>%
  #removing treatments that were low levels of resource manipulation, not manipulating any resource, pulses
  #plus filtering to get only treatments
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")

divCDRe001<-div%>%
  filter(site_code=="CDR"&treatment==1|treatment==6|treatment==8|treatment==9,plot_mani>0)
divCDRe002<-div%>%
  filter(site_code=="CDR"&treatment=='1_f_u_n'|treatment=='6_f_u_n'|treatment=='8_f_u_n'|treatment=='9_f_u_n',plot_mani>0)
  

divTrt<-rbind(divTrt1, divCDRe002, divCDRe001)

##16% of our data is from CDR and 10% is from KNZ
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
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))

#full dataset
write.csv(ForAnalysis, "ForBayesianAnalysis_May2017.csv")

#8 yr or less
ForAnalysis8yr <- ForAnalysis%>%
  filter(treatment_year<9)

##18% of our our data is from CDR and KNZ

write.csv(ForAnalysis8yr, "ForBayesianAnalysis_8yr_May2017.csv")


#Plot of 9 year data used in paper.
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



#8 year, single resource manipulations only
resource8 <- ForAnalysis8yr%>%
  left_join(expInfo)%>%
  #filter out anything with more than one resource manipulated
  filter(plot_mani<2)%>%
  #make resources binary
  mutate(n1=ifelse(n>0, 1, 0), p1=ifelse(p>0, 1, 0), other_nut=ifelse(k>0, 1, ifelse(other_trt=='mirconutrients and lime added', 1, ifelse(other_trt=='lime added', 1, 0))), CO2_1=ifelse(CO2>0, 1, 0), irr=ifelse(precip>1, 1, 0), drought=ifelse(precip<0, 1, 0), sum=n1+p1+other_nut+CO2_1+irr+drought)%>%
  #drop megarich because they add NPK to all megaliths, so no true resource control
  filter(project_name!='MEGARICH', sum==1)%>%
  mutate(resource=ifelse(n1==1, 'n', ifelse(p1==1, 'p', ifelse(other_nut==1, 'other_nut', ifelse(CO2_1==1, 'CO2', ifelse(irr==1, 'irrigation', 'drought'))))))

write.csv(resource8, 'ForBayesianAnalysis_singleresource_May2017.csv')
  
  
  


#9+ year datasets (all years)
ForAnalysis9yr <- ForAnalysis%>%
  filter(experiment_length>8)
write.csv(ForAnalysis9yr, "ForBayesianAnalysis_9plusyr_May2017.csv")

# #absolute value
# ForAnalysisAbsValue <- ForAnalysis9yr%>%
#   mutate(mean_change=abs(mean_change), dispersion_change=abs(dispersion_change), H_change=abs(H_change), SimpEven_change=abs(SimpEven_change), S_PC=abs(S_PC))
# write.csv(ForAnalysisAbsValue, "ForBayesianAnalysis_abs value_9yr_Dec2016.csv")




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
<<<<<<< HEAD
=======
  write.csv(ForANPPAnalysis, "ForBayesianAnalysisANPP_Oct2017.csv")
  
ForANPPAnalysis9yr<-ForANPPAnalysis%>%
  filter(treatment_year<10)
write.csv(ForANPPAnalysis9yr, "ForBayesianAnalysisANPP_9yr_Dec2016.csv")
>>>>>>> 455efdf31caabc67e5a48e71fa15195d60f20a05

write.csv(ForANPPAnalysis, "ForBayesianAnalysisANPP_Oct2017.csv")

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
