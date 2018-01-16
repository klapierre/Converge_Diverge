library(tidyverse)

setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

expInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  select(-X, -public)
siteInfo <- read.csv('SiteExperimentDetails_Dec2016.csv')%>%
  select(-X, -public)
latLong <- read.csv('siteList_LatLong.csv')%>%
  mutate(site_code=name)%>%
  select(-name)

adamMetadata <- expInfo%>%
  left_join(siteInfo)%>%
  left_join(latLong)%>%
  filter(project_name=='watering'|project_name=='MAT2'|project_name=='MNT'|project_name=='PQ'|project_name=='BioCON'|project_name=='e001'|project_name=='e002'|project_name=='GCE'|project_name=='T7'|project_name=='BFFert'|project_name=='IRG'|project_name=='pplots'|project_name=='E6'|project_name=='246Nfert'|community_type=='DryBowman'|project_name=='CXN'|site_code=='TAS')%>%
  select(-calendar_year, -treatment_year, -pulse, -resource_mani, -plot_mani, -max_trt, -factorial, -id, -community_type)%>%
  unique()

adamMetadata2 <- adamMetadata%>%
  select(site_code, project_name, MAP, MAT, rrich, experiment_length, anpp, latitude, longitude)%>%
  unique()

#write adamMetadata2 and add trt column by hand
# write.csv(adamMetadata2, 'Langley_metadata_01162017.csv')