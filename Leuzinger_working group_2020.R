library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

expInfo <- read.csv('ExperimentInformation_March2019.csv')%>%
  select(-X)

siteData <- read.csv('SiteExperimentDetails_March2019.csv')%>%
  select(-X)%>%
  rename(site_productivity=anpp)

anppPublic <- read.csv('ANPP_Oct2017.csv')%>%
  select(-X)%>%
  left_join(siteData)%>%
  left_join(expInfo)%>%
  filter(public==1)

write.csv(anppPublic, 'CoRRE_anpp_public_Feb2020.csv', row.names=F)
