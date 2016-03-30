#kim's
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#meghan's
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)

#read in site data
siteExperimentDetails <- read.csv('SiteExperimentDetails_Feb2016a.csv')%>%
  select(site_code, project_name, community_type, public)

#read in species data
sppAbund <- read.csv('SpeciesRawAbundance_Feb2016.csv')%>%
  select(-X)%>%
  left_join(siteExperimentDetails)%>%
  filter(public==1)

#write species data
write.csv(sppAbund, 'SpeciesRawAbundance_Feb2016_public.csv')



#read in ANPP data
ANPP <- read.csv('ANPP_Feb2016.csv')%>%
  select(-X)%>%
  left_join(siteExperimentDetails)%>%
  filter(public==1)

#write ANPP data
write.csv(ANPP, 'ANPP_Feb2016_public.csv')
