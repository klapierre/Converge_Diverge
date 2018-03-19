library(lme4)
library(lsmeans)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

###read in data

#experiment information
expInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep='_'))

anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)

#diversity data
div <- read.csv('DiversityMetrics_Nov2017.csv')%>%
  separate(exp_year, c('site_code', 'project_name','community_type', 'calendar_year'), sep='::')%>%
  mutate(calendar_year=as.integer(calendar_year))%>%
  left_join(expInfo)%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  #create e^H metric
  mutate(expH=exp(H))

#anpp data
anpp<-read.csv("ANPP_Oct2017.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)

#appearance/disappearance data
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

###calculate change in dispersion, H, S, and evenness
#subset out controls and treatments
divControls <- subset(div, subset=(plot_mani==0))%>%
  select(exp_year, dispersion, expH, S, SimpEven, calendar_year, treatment_year)
  names(divControls)[names(divControls)=='dispersion'] <- 'ctl_dispersion'
  names(divControls)[names(divControls)=='expH'] <- 'ctl_expH'
  names(divControls)[names(divControls)=='S'] <- 'ctl_S'
  names(divControls)[names(divControls)=='SimpEven'] <- 'ctl_SimpEven'
divTrt1 <- div%>%
  #filtering to get only non-e002, e001 treatments
  filter(pulse==0, plot_mani>0, project_name!="e001"&project_name!="e002")

#calculate average dispersion among just control plots (this is an estimate of average community dissimilarity, to use as a baseline for dissimilarity change between treatment and control plots)
summary(divControls$ctl_dispersion) #mean=0.2811, median=0.2778


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
#calculate change in disperion, expH, S, and evenness
  mutate(dispersion_change=dispersion-ctl_dispersion, 
         expH_PC=(expH-ctl_expH)/ctl_expH, 
         S_PC=(S-ctl_S)/ctl_S, 
         SimpEven_change=SimpEven-ctl_SimpEven)%>%
  select(exp_year, treatment_year, treatment, plot_mani, mean_change, dispersion_change, expH_PC,  SimpEven_change, S_PC, site_code, project_name, community_type, calendar_year)

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


###merging with experiment (treatment) information
SiteExp<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  select(-X)

ForAnalysis<-merge(divCompare, SiteExp, by=c("site_code","project_name","community_type"))

###read in pairwise multivariate distance
pair <- read.csv('Bray_Curtis_Ave_dissim_03162018.csv')%>%
  select(-X)%>%
  left_join(expInfo)%>%
  filter(plot_mani==0)%>%
  select(site_code, project_name, community_type, treatment2, calendar_year, BC_between_diff)%>%
  mutate(treatment=treatment2)%>%
  select(-treatment2)

ForAnalysis <- ForAnalysis%>%
  left_join(ForAnalysis)%>%
  filter(!is.na(BC_between_diff))


#full dataset
# write.csv(ForAnalysis, "ForBayesianAnalysis_May2017.csv")

#8 yr or less
ForAnalysis8yr <- ForAnalysis%>%
  filter(treatment_year<9)
#subset out datasets with less than 5 temporal data points
numPoints <- ForAnalysis8yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
ForAnalysis8yrPoints <- ForAnalysis%>%
  left_join(numPoints)%>%
  filter(num_datapoints>2)

##18% of our our data is from CDR and KNZ

# write.csv(ForAnalysis8yrPoints, "ForBayesianAnalysis_8yr_May2017.csv")


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

#9+ year datasets (all years)
ForAnalysis9yr <- ForAnalysis%>%
  filter(experiment_length>8)
# write.csv(ForAnalysis9yr, "ForBayesianAnalysis_9plusyr_May2017.csv")

#absolute value
ForAnalysisAbsValue <- ForAnalysis9yr%>%
  mutate(mean_change=abs(mean_change), dispersion_change=abs(dispersion_change), expH_PC=abs(expH_PC), SimpEven_change=abs(SimpEven_change), S_PC=abs(S_PC))
# write.csv(ForAnalysisAbsValue, "ForBayesianAnalysis_abs value_9yr_Dec2016.csv")


# ###getting the treatment interaction types -- used for past analyses and ANPP analysis
# trtType <- ForAnalysis%>%
#   left_join(expInfo)%>%
#   #create drought and irrigation categories
#   mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))%>%
#   #create categorical treatment type column
#   mutate(trt_type=ifelse(plot_mani==2&n>0&drought<0, 'N+drought', ifelse(plot_mani==2&n>0&irrigation>0, 'N+irr', ifelse(plot_mani==2&n>0&p>0, 'N+P', ifelse(plot_mani==2&p>0&k>0, 'P+K', ifelse(plot_mani==2&n>0&CO2>0, 'N+CO2', ifelse(plot_mani==2&CO2>0&irrigation>0, 'CO2+irr', ifelse(plot_mani==3&n>0&p>0&k>0, 'N+P+K', ifelse(plot_mani==3&n>0&CO2>0&irrigation>0, 'N+CO2+irr', ifelse(plot_mani==4&n>0&p>0&k>0&irrigation>0, 'N+P+K+irr', ifelse(plot_mani==1&n>0, 'N', ifelse(plot_mani==1&p>0, 'P', ifelse(plot_mani==1&irrigation>0, 'irr', ifelse(plot_mani==1&drought<0, 'drought', ifelse(plot_mani==1&CO2>0, 'CO2', ifelse(n==0&drought==0&irrigation==0&p==0&k==0&CO2==0, 'other', 'resource+other'))))))))))))))))%>%
#   #drop experiments that we can't run the analyses on
#   filter(plot_mani<6, treatment_year!=0, anpp!='NA')%>%
#   #keep just relevent column names for this analysis
#   select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, dispersion_change, SimpEven_change, S_PC, experiment_length, resource_other)%>%
#   #keep final year only
#   group_by(site_code, project_name, community_type, treatment)%>%
#   filter(treatment_year==max(treatment_year))%>%
#   ungroup()

# write.csv(trtType, 'treatment interactions_11152017.csv')

###getting resource*non-resource interactions
###4 analyses: (1) single resource, (2) single non-resource, (3) 2-way interactions, (4) 3+ way interactions
#analysis 1: single resource
singleResource <- ForAnalysis%>%
  select(-plot_mani)%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==1, plot_mani==1)%>%
  #set CEH Megarich nutrient values to 0 (added to all megaliths, not a treatment)
  mutate(n2=ifelse(site_code=='CEH', 0, n), p2=ifelse(site_code=='CEH', 0, p), k2=ifelse(site_code=='CEH', 0, k))%>%
  #drop lime added, as only one trt does this
  filter(other_trt!='lime added')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(n2>0, 'N', ifelse(p2>0, 'P', ifelse(k2>0, 'K', ifelse(precip<0, 'drought', ifelse(precip>0, 'irr', ifelse(CO2>0, 'CO2', 'precip_vari')))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, BC_between_diff, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

# singleResource8yr <- singleResource%>%
#   filter(treatment_year<9)
# # write.csv(singleResource8yr, 'ForAnalysis_singleResource8yr.csv')
# singleResourceAbs <- singleResource%>%
#   mutate(S_PC_abv=abs(S_PC))%>%
#   select(-S_PC)
# # write.csv(singleResourceAbs, 'ForAnalysis_singleResourceAbs.csv')
# singleResource9yr <- singleResource%>%
#   filter(experiment_length>8)
# # write.csv(singleResource9yr, 'ForAnalysis_singleResource9yr.csv')

# ##check the treatment designations are correct; only works if you don't run select line in previous step
# temp <- singleResource%>%
#   select(site_code, project_name, treatment, trt_type, n, p, k, CO2, precip, other_trt, mow_clip, burn, herb_removal, temp, plant_mani, plot_mani)%>%
#   group_by(site_code, project_name, treatment, trt_type)%>%
#   unique()%>%
#   ungroup()

#analysis 2: single non-resource
singleNonresource <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(resource_mani==0, plot_mani==1)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==1, 'burn', ifelse(mow_clip==1, 'mow_clip', ifelse(herb_removal==1, 'herb_rem', ifelse(temp>0, 'temp', ifelse(plant_trt==1, 'plant_mani', 'other'))))))%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, BC_between_diff, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

# singleNonresource8yr <- singleNonresource%>%
#   filter(treatment_year<9)
# # write.csv(singleNonresource8yr, 'ForAnalysis_singleNonresource8yr.csv')
# singleNonresourceAbs <- singleNonresource%>%
#   mutate(S_PC_abv=abs(S_PC))%>%
#   select(-S_PC)
# # write.csv(singleNonresourceAbs, 'ForAnalysis_singleNonresourceAbs.csv')
# singleNonresource9yr <- singleNonresource%>%
#   filter(experiment_length>8)
# # write.csv(singleNonresource9yr, 'ForAnalysis_singleNonresource9yr.csv')

###check the treatment designations are correct; only works if you don't run select line in previous step
# temp <- singleNonresource%>%
#   select(site_code, project_name, treatment, trt_type, n, p, k, CO2, precip, other_trt, mow_clip, burn, herb_removal, temp, plant_mani, plot_mani)%>%
#   group_by(site_code, project_name, treatment, trt_type)%>%
#   unique()%>%
#   ungroup()
  
#analysis 3: 2-way interactions
twoWay <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani==2)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(resource_mani==1&burn==1, 'R*burn', ifelse(resource_mani==1&mow_clip==1, 'R*mow_clip', ifelse(resource_mani==1&herb_removal==1, 'R*herb_rem', ifelse(resource_mani==1&temp>0, 'R*temp', ifelse(resource_mani==1&plant_trt==1, 'R*plant_mani', ifelse(resource_mani==1&other_trt!=0, 'R*other', ifelse(n>0&p>0, 'R*R', ifelse(n>0&CO2>0, 'R*R', ifelse(n>0&precip!=0, 'R*R', ifelse(p>0&k>0, 'R*R', ifelse(CO2>0&precip!=0, 'R*R', 'N*N'))))))))))))%>%
  #drop R*herb_removal (single rep)
  filter(trt_type!='R*herb_rem')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, BC_between_diff, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

# twoWay8yr <- twoWay%>%
#   filter(treatment_year<9)
# # write.csv(twoWay8yr, 'ForAnalysis_twoWay8yr.csv')
# twoWayAbs <- twoWay%>%
#   mutate(S_PC_abv=abs(S_PC))%>%
#   select(-S_PC)
# # write.csv(twoWayAbs, 'ForAnalysis_twoWayAbs.csv')
# twoWay9yr <- twoWay%>%
#   filter(experiment_length>8)
# # write.csv(twoWay9yr, 'ForAnalysis_twoWay9yr.csv')

# ##check the treatment designations are correct; only works if you don't run select line in previous step
# temp <- twoWay%>%
#   select(site_code, project_name, treatment, trt_type, n, p, k, CO2, precip, other_trt, mow_clip, burn, herb_removal, temp, plant_mani, plot_mani, resource_mani)%>%
#   group_by(site_code, project_name, treatment, trt_type)%>%
#   unique()%>%
#   ungroup()


#analysis 4: 3+ way interactions
threeWay <- ForAnalysis%>%
  left_join(expInfo)%>%
  #filter pretrt data
  filter(treatment_year!=0)%>%
  #filter just single resource manipulations
  filter(plot_mani>2, plot_mani<6)%>%
  #drop tilled, stone, and fungicide as only one trt each do these; drop 'current pattern' of rainfall, as this is basically a control
  filter(other_trt!='tilled', other_trt!='shallow soil', other_trt!='fungicide added', other_trt!='current pattern')%>%
  #create categorical treatment type column
  mutate(trt_type=ifelse(burn==0&mow_clip==0&herb_removal==0&temp==0&plant_trt==0, 'all_resource', ifelse(n==0&p==0&k==0&CO2==0&precip==0, 'all_nonresource', 'both')))%>%
  #drop single all-nonresource treatment (NIN herbdiv 5NF)
  filter(trt_type!='all_nonresource')%>%
  #keep just relevent column names for this analysis
  select(site_code, project_name, community_type, exp_year, treatment_year, calendar_year, treatment, trt_type, mean_change, BC_between_diff, S_PC, expH_PC, experiment_length, rrich, anpp, MAT, MAP)

# threeWay8yr <- threeWay%>%
#   filter(treatment_year<9)
# # write.csv(threeWay8yr, 'ForAnalysis_threeWay8yr.csv')
# threeWayAbs <- threeWay%>%
#   mutate(S_PC_abv=abs(S_PC))%>%
#   select(-S_PC)
# # write.csv(threeWayAbs, 'ForAnalysis_threeWayAbs.csv')
# threeWay9yr <- threeWay%>%
#   filter(experiment_length>8)
# # write.csv(threeWay9yr, 'ForAnalysis_threeWay9yr.csv')

###check the treatment designations are correct; only works if you don't run select line in previous step
# temp <- threeWay%>%
#   select(site_code, project_name, treatment, trt_type, n, p, k, CO2, precip, other_trt, mow_clip, burn, herb_removal, temp, plant_mani, plot_mani, resource_mani)%>%
#   group_by(site_code, project_name, treatment, trt_type)%>%
#   unique()%>%
#   ungroup()


#combine for analysis - one big model, 19 trt types
allAnalysis <- rbind(singleResource, singleNonresource, twoWay, threeWay)
# write.csv(allAnalysis, 'ForAnalysis_allAnalysis.csv')
# allAnalysis8yr <- allAnalysis%>%
#   filter(treatment_year<9)
# # write.csv(allAnalysis8yr, 'ForAnalysis_allAnalysis8yr.csv')
# allAnalysisAbs <- allAnalysis8yr%>%
#   mutate(S_PC_abv=abs(S_PC))%>%
#   select(-S_PC)
# # write.csv(allAnalysisAbs, 'ForAnalysis_allAnalysisAbs.csv')
# allAnalysis9yr <- allAnalysis%>%
#   filter(experiment_length>8)
# # write.csv(allAnalysis9yr, 'ForAnalysis_allAnalysis9yr.csv')

#subset out datasets with less than 3 temporal data points
numPoints <- allAnalysis%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
#which are dropped? GVN FACE
allAnalysisAllDatasets <- allAnalysis%>%
  left_join(numPoints)%>%
  filter(num_datapoints>2)
# write.csv(allAnalysisAllDatasets, 'ForAnalysis_allAnalysisAllDatasets.csv')


#subset out datasets 10 years and shorter
allAnalysis10yr <- allAnalysisAllDatasets%>%
  filter(treatment_year<11)%>%
  select(-num_datapoints)
numPoints <- allAnalysis10yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
#which are dropped? GVN FACE, KNZ BGP, SEV Nfert
allAnalysis10yrPoints <- allAnalysis10yr%>%
  left_join(numPoints)%>%
  filter(num_datapoints>2)
# write.csv(allAnalysis10yrPoints, 'ForAnalysis_allAnalysis10yr.csv')

#subset out datasets 15 years and shorter
allAnalysis15yr <- allAnalysisAllDatasets%>%
  filter(treatment_year<16)%>%
  select(-num_datapoints)
numPoints <- allAnalysis15yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
#nothing would be dropped
# write.csv(allAnalysis15yr, 'ForAnalysis_allAnalysis15yr.csv')

#subset out datasets 20 years and shorter
allAnalysis20yr <- allAnalysisAllDatasets%>%
  filter(treatment_year<21)%>%
  select(-num_datapoints)
numPoints <- allAnalysis20yr%>%
  select(site_code, project_name, community_type, treatment, treatment_year)%>%
  unique()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(num_datapoints=length(treatment_year))
#nothing would be dropped
# write.csv(allAnalysis20yr, 'ForAnalysis_allAnalysis20yr_pairwise.csv')

#subset out 20th or final year of all data
allAnalysisFinalYear <- allAnalysis20yr%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  filter(treatment_year==max(treatment_year))
# write.csv(allAnalysisFinalYear, 'ForAnalysis_allAnalysisFinalYear.csv')


#get magnitudes
allAnalysisMag <- allAnalysis20yr%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  filter(treatment_year==max(treatment_year))%>%
  left_join(expInfo)

Nmag <- allAnalysisMag%>%
  filter(project_name!='MEGARICH')%>%
  filter(n>0)
# write.csv(Nmag, 'ForAnalysis_allAnalysisNmag.csv')
H2OMag <- allAnalysisMag%>%
  filter(precip!=0)
# write.csv(H2OMag, 'ForAnalysis_allAnalysisH2Omag.csv')



#a few prelim figures
allAnalysis8yrFig <- allAnalysis10yrPoints%>%
  mutate(temp=paste(site_code, project_name, community_type, treatment))%>%
  filter(treatment_year<8, experiment_length>8)
eight <- ggplot(data=allAnalysis8yrFig, aes(x=treatment_year, y=mean_change, color=temp)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F) +
  geom_point() +
  theme(legend.position='none') +
  facet_wrap(~temp)


allAnalysis10yrFig <- allAnalysis10yrPoints%>%
  mutate(temp=paste(site_code, project_name, community_type, treatment))%>%
  filter(experiment_length>8)
ten <- ggplot(data=allAnalysis10yrFig, aes(x=treatment_year, y=mean_change, color=temp)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F) +
  geom_point() +
  theme(legend.position='none') +
  facet_wrap(~temp)

grid.arrange(eight, ten, ncol=2)

allAnalysis15yrFig <- allAnalysis15yr%>%
  mutate(temp=paste(site_code, project_name, community_type, treatment))
ggplot(data=allAnalysis15yrFig, aes(x=treatment_year, y=mean_change, color=temp)) +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F) +
  theme(legend.position='none')

allAnalysis20yrFig <- allAnalysis20yr%>%
  mutate(temp=paste(site_code, project_name, community_type, treatment))%>%
  filter(experiment_length>15)
ggplot(data=allAnalysis20yr, aes(x=treatment_year, y=BC_between_diff, color=trt_type)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F) +
  theme(legend.position='none') +
  facet_wrap(~trt_type)
ggplot(data=allAnalysis20yr, aes(x=treatment_year, y=mean_change, color=trt_type)) +
  geom_point() +
  geom_smooth(method='lm', formula=y~x+I(x^2), se=F) +
  theme(legend.position='none') +
  facet_wrap(~trt_type)

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
