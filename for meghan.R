library(dplyr)
library(tidyr)
library(plyr)

setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\to yang')

test <- read.csv('CORRE_raw_abundance.csv')
trt <- read.csv('CORRE_treatment_summary.csv')%>%
  select(-X)

cdrOnly <- test%>%
  subset(project_name=='e001'&community_type=='D')

irrOnly <- test%>%
  subset(project_name=='IRG'&community_type=='l')

bothOnly <- cdrOnly%>%
  rbind(irrOnly)%>%
  select(-X)%>%
  left_join(trt)%>%
  select(-nutrients, -light, -water, -carbon, -other_manipulation, -CO2, -temp, -mow_clip, -burn, -herb_removal, -trt_details, -successional, -plant_mani, -pulse, -resource_mani, -max_trt, -public)

write.csv(bothOnly, 'CORRE_irrigation_e001D_subset.csv')



jasperOnly <- test%>%
  filter(site_code=='JSP')%>%
  left_join(trt)%>%
  select(-nutrients, -light, -water, -carbon, -other_manipulation, -p, -k, -mow_clip, -burn, -herb_removal, -trt_details, -other_trt, -successional, -plant_mani, -pulse, -resource_mani, -public)


write.csv(jasperOnly, 'CORRE_jasper_subset.csv')
