library(tidyverse)

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#meghan's
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

dat<-read.csv("SpeciesRawAbundance_Oct2017.csv")%>%
  group_by(site_code, project_name, community_type)%>%
  mutate(maxyear=max(treatment_year))%>%
  filter(maxyear>9)%>%
  select(-X)

trt<- read.csv('ExperimentInformation_May2017.csv')%>%
  select(-X)

dat2<-merge(trt, dat, by=c("site_code","project_name","community_type","treatment", "calendar_year","treatment_year"))%>%
  select(-nutrients, -light, -carbon, -water, -other_manipulation, -max_trt, -public, -factorial, -block, -maxyear)

dat_touse<-dat2%>%
  filter(plot_mani!=0)%>%
  group_by(site_code, project_name, community_type, plot_mani)%>%
  mutate(mintrt=min(plot_mani))%>%
  filter(mintrt==1, plot_mani==1|plot_mani==0)%>%
  ungroup()%>%
  select(site_code, project_name, community_type)%>%
  unique()

datsub<-merge(dat2, dat_touse, by=c("site_code","project_name",'community_type'))%>%
  filter(plot_mani==0|plot_mani==1)

#get the arc mnt, mat2, cdr e001, e002 to paste into old dataset
mnt<-dat2%>%
  filter(project_name=='MNT')

mat2<-dat2%>%
  filter(project_name=='MAT2')

e001<-dat2%>%
  filter(project_name=='e001')%>%
  filter(treatment==1|treatment==5|treatment==7|treatment==9)

e002<-dat2%>%
  filter(project_name=='e002')%>%
  filter(treatment=="1_f_u_n"|treatment=="5_f_u_n"|treatment=="7_f_u_n"|treatment=="9_f_u_n")

extras <- rbind(mnt, mat2, e001, e002)

foradam<-rbind(datsub, extras)

write.csv(foradam, 'CORRE_tenyears_extras_raw.csv')
