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

rawData <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_spatial\\anpp_spatial_data.csv')# get SD and means and use this backtransform the chains

mean(rawData$anpp_sp_cv, na.rm=T) #29.78862
sd(rawData$anpp_sp_cv, na.rm=T) #16.73752


#treatment info
trtInfo<-read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_spatial\\site_list.csv')

################################################################################
################################################################################
#only run to generate initial chains files
#raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_spatial\\anpp_spatial_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_spatial\\anpp_spatial_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_spatial\\anpp_spatial_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_spatial\\anpp_spatial_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsANPP <- rbind(chains1, chains2, chains3, chains4)
# 
# write.csv(chainsANPP, "ANPP_Spatial_Chains.csv")
chainsANPP<-read.csv("ANPP_Spatial_Chains.csv")
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsANPP$mu.1))
# plot(density(chainsANPP$mu.2))
# plot(density(chainsANPP$mu.3))
# 
# 
# #get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
# chainsANPP2 <- chainsANPP%>%
#   select(lp__,
#          #plot_mani intercepts (center digit): 1=plot mani 1, 2=plot mani 2, 3=plot mani 3, 4=plot mani 4, 5= plot mani 5
#          U.1.1,
#          U.2.1,
#          U.3.1,
#          U.4.1,
#          U.5.1,
#          #plot_mani intercepts (center digit): 1=plot mani 1, 2=plot mani 2, 3=plot mani 3, 4=plot mani 4, 5= plot mani 5
#          U.1.2,
#          U.2.2,
#          U.3.2,
#          U.4.2,
#          U.5.2,
#          #plot_mani intercepts (center digit): 1=plot mani 1, 2=plot mani 2, 3=plot mani 3, 4=plot mani 4, 5= plot mani 5
#          U.1.3,
#          U.2.3,
#          U.3.3,
#          U.4.3,
#          U.5.3,
#          #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
#          D.1.1,
#          D.1.2,
#          D.1.3,
#          #richness intercept, linear, and quad slopes (center digit): 2=richness
#          D.2.1,
#          D.2.2,
#          D.2.3,
#          #MAP intercept, linear, and quad slopes (center digit): 1=MAP
#          E.1.1,
#          E.1.2,
#          E.1.3,
#          #MAT intercept, linear, and quad slopes (center digit): 2=MAT
#          E.2.1,
#          E.2.2,
#          E.2.3,
#          #overall intercept, linear, and quad slopes
#          mu.1, mu.2, mu.3)%>%
#   gather(key=parameter, value=value, U.1.1:mu.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median2=ifelse(diff==-2, 0, median))

write.csv(chainsANPP2, 'bayesian_ANPP_spatial_output_summary_10212017.csv')

chainsANPP2 <- read.csv('bayesian_ANPP_spatial_output_summary_10212017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsANPP[,2487:3551]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsANPP[, 2487:3551])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsANPP[,2487:3551]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsANPP[, 2487:3551])'] <- 'sd'
# 
# chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
#   #split names into parts
#   separate(parameter, c('B', 'id', 'parameter'))%>%
#   select(-B)%>%
#   #rename parts to be more clear
#   mutate(parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          id=as.integer(id))%>%
#   #if 95% confidence interval overlaps 0, then set mean to 0
#   mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
#   #spread by variable
#   select(id, parameter, mean)%>%
#   spread(key=parameter, value=mean)
# 
# write.csv(chainsFinal, 'bayesian_anpp_spatial_output_mean sd_102117.csv')

chainsFinal <- read.csv('bayesian_anpp_spatial_output_mean sd_102117.csv')

#merge together with experiment list
chainsFinal <- chainsFinal%>%
  arrange(id)

chainsExperiment<-cbind(chainsFinal, trtInfo)

minyear<-chainsExperiment%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(min_year=min(treatment_year))


mean<-chainsExperiment%>%
  left_join(minyear)%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  #get estimates at 10 years
  mutate(final_year_estimate=(intercept + alt_length*linear + (alt_length^2)*quadratic)*0.6112685 + 0.362312)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,',
         curve5='), colour=gray) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, experiment_length, curve5, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

write.csv(mean, "anpp_to_clean.csv", row.names = F)

#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  #coord_cartesian(xlim=c(0,24), ylim=c(-2,4))  +
  #scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  #ylim(-3,6) +
  xlab('Standardized Year') +
  ylab('ANPP Percent Change')

meanPlot <- meanPlot + 
#below are the individual treatment lines
stat_function(fun=function(x){}, size=3, xlim=c(0,23), colour='black')
print(meanPlot)
#export at 1200x1000


###by resource mani




###summary stats from bayesian output
chainsANPP3<-chainsANPP2%>%
  separate(parameter, c('Type', 'id', 'parama'),remove=F)%>%
 #rename parts to be more clear
 mutate(parameter2=ifelse(parama==1, 'intercept', ifelse(parama==2, 'linear', 'quadratic')))%>%
  mutate(param=paste(Type, id, sep="."))%>%
  mutate(CI=sd*2)

int<-chainsANPP3%>%
  filter(parameter2=="intercept"|parameter=="mu.1")
lin<-chainsANPP3%>%
  filter(parameter2=="linear"|parameter=="mu.2")
quad<-chainsANPP3%>%
  filter(parameter2=="quadratic"|parameter=="mu.3")

#mean plots
ggplot(data=int, aes(x=param, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('E.1', 'E.2', 'D.1', 'D.2', 'U.5', 'U.4', 'U.3', 'U.2','U.1','mu.1'), labels=c('MAP','MAT','ANPP','Gamma Diversity','5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations','1 Manipulation','Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

ggplot(data=lin, aes(x=param, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('E.1', 'E.2', 'D.1', 'D.2', 'U.5', 'U.4', 'U.3', 'U.2','U.1','mu.2'), labels=c('MAP','MAT','ANPP','Gamma Diversity','5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations','1 Manipulation','Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  ylim(-.15, .15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

ggplot(data=quad, aes(x=param, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('E.1', 'E.2', 'D.1', 'D.2', 'U.5', 'U.4', 'U.3', 'U.2','U.1','mu.3'), labels=c('MAP','MAT','ANPP','Gamma Diversity','5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations','1 Manipulation','Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  ylim(-.015, .015) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()
