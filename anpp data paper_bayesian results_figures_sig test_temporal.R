library(ggplot2)
library(grid)
library(mgcv)
library(plyr)
library(dplyr)
library(tidyr)

#kim's laptop
setwd('C:\\Users\\Kim\\Desktop\\bayesian output')

#meghan's desktop
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

##################################################################################
##################################################################################
#experiment information
expInfo <- read.csv('ExperimentInformation_ANPP_Oct2017.csv')

expInfo2 <- read.csv('ExperimentInformation_May2017.csv')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), precip=mean(precip))

rawData <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_temporal\\anpp_temporal_data.csv')# get SD and means and use this backtransform the chains

mean(rawData$anpp_temp_cv, na.rm=T) #38.2354
sd(rawData$anpp_temp_cv, na.rm=T) #15.84088


#treatment info THIS FILES DOES NOT EXIST
trtInfo<-read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_temporal\\site_list.csv')

###############################################################################
###############################################################################
#only run to generate initial chains files
#raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_temporal\\anpp_temporal_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_temporal\\anpp_temporal_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_temporal\\anpp_temporal_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_temporal\\anpp_temporal_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsANPP <- rbind(chains1, chains2, chains3, chains4)
# 
# write.csv(chainsANPP, "ANPP_temporal_chains.csv")
chainsANPP<-read.csv("ANPP_temporal_chains.csv")

#density plot of chains --------------------------------------------------------
plot(density(chainsANPP$mu))#I would say this only looks okay



#get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
chainsANPP2 <- chainsANPP%>%
  select(lp__,
         #plot_mani output (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1,
         U.2,
         U.3,
         U.4,
         U.5,
        #ANPP output,
         D.1,
         #richness output, l
         D.2,
         #MAP output,
         E.1,
         #MAT intercept,
         E.2,
         #overall intercept, linear, and quad slopes
         mu)%>%
  gather(key=parameter, value=value, U.1:mu)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value), lower=quantile(value, 0.025), upper=quantile(value, 0.975))%>%
  mutate(lower2=median-2*sd, upper2=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median2=ifelse(diff==-2, 0, median))

write.csv(chainsANPP2, 'bayesian_ANPP_temporal_output_summary_10212017.csv')

chainsANPP2 <- read.csv('bayesian_ANPP_temporal_output_summary_10212017.csv')

###summary stats from bayesian output
chainsANPP2<-chainsANPP2%>%
  mutate(CI=sd*2)

theme_set(theme_bw(12))

#stat plot
ggplot(data=chainsANPP2, aes(x=parameter, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('E.1', 'E.2', 'D.1', 'D.2', 'U.5', 'U.4', 'U.3', 'U.2','U.1','mu'), labels=c('MAP','MAT','ANPP','Gamma Diversity','5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations','1 Manipulation','Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()
