library(tidyr)
library(dplyr)

#kim's
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')


###read in data

#experiment information
expInfo <- read.csv('ExperimentInformation_Nov2015.csv')%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep='::'))

#diversity data
div <- read.csv('DiversityMetrics_Nov2015.csv')



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














