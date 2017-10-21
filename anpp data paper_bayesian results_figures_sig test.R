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

rawData <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_mean_change\\anpp_mean_change_data.csv')# get SD and means and use this backtransform the chains

mean(rawData$anpp_PC, na.rm=T) #0.3760375
sd(rawData$anpp_PC, na.rm=T) #0.6968239


#treatment info
trtInfo<-read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_mean_change\\site_list.csv')

################################################################################
################################################################################
# #only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_mean_change\\anpp_mean_change_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_mean_change\\anpp_mean_change_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_mean_change\\anpp_mean_change_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\nate_results\\anpp\\anpp_mean_change\\anpp_mean_change_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsANPP <- rbind(chains1, chains2, chains3, chains4)
# 
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
#          #plot_mani intercepts (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.1,
#          U.2.1,
#          U.3.1,
#          U.4.1,
#          #plot_mani linear slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.2,
#          U.2.2,
#          U.3.2,
#          U.4.2,
#          #plot_mani quad slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.3, 
#          U.2.3, 
#          U.3.3,
#          U.4.3, 
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
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))
# 
# write.csv(chainsANPP2, 'bayesian_ANPP_output_summary_10162017.csv')

chainsANPP2 <- read.csv('bayesian_ANPP_output_summary_10162017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsANPP[,1866:2621]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsANPP[, 1866:2621])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsANPP[,1866:2621]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsANPP[, 1866:2621])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_ampp_output_mean sd_101617.csv')

chainsFinal <- read.csv('bayesian_ampp_output_mean sd_101617.csv')

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
stat_function(fun=function(x){(0.7048927463627 + 0*x + 0.026608528145043*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,13), colour="gray") +
  stat_function(fun=function(x){(0 + -0.293589368291105*x + 0.0287403453996*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,13), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0.596602024874575 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-0.63139246705765 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.641276376937955 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.6226030607248 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.62162233424415 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.76015482935905 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.7498961622133 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.79566742154 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.703938635991 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.71521147477085 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.6876313164489 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0.78486963985015 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0.6052342217654 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,31), colour="gray") +
  stat_function(fun=function(x){(1.3300086605 + -0.18655944646*x + 0.008209384191*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,31), colour="gray") +
  stat_function(fun=function(x){(1.66614955355 + -0.18666291198*x + 0.0080198254585*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,31), colour="gray") +
  stat_function(fun=function(x){(-1.099726676755 + 0.1016195028034*x + -0.0033692678014685*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,31), colour="gray") +
  stat_function(fun=function(x){(0.713251491412205 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,31), colour="gray") +
  stat_function(fun=function(x){(1.11984465655 + 0.098904062271895*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,31), colour="gray") +
  stat_function(fun=function(x){(-0.744864790783305 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,30), colour="gray") +
  stat_function(fun=function(x){(0.925188235165 + 0*x + 0.0029828005394151*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,30), colour="gray") +
  stat_function(fun=function(x){(1.18647473555 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,30), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,27), colour="gray") +
  stat_function(fun=function(x){(1.62650924845 + -0.1599998247715*x + 0.003751217809125*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,27), colour="gray") +
  stat_function(fun=function(x){(2.384622966 + -0.10037026209075*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,27), colour="gray") +
  stat_function(fun=function(x){(0 + -0.531441939975*x + 0.0483915845951*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(2.06003258645 + -0.6875259241*x + 0.0540857237885*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(2.319305386 + -0.6793659639*x + 0.053559244375*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0.76526394340594 + -0.83685356325*x + 0.08270921852*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(2.508355223 + -0.97691600015*x + 0.088232157285*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(2.7133112305 + -0.98495289475*x + 0.088063644315*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + -0.3147085768448*x + 0.035370739611705*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(1.3734693441 + -0.47076246739*x + 0.0412748756519465*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(1.77850418745 + -0.4303190650245*x + 0.04001803513225*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0.656730841225345 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(-0.671287825878318 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.570909330482553 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.6321712379237 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.7667240208568 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.7102781741615 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.75273408702125 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.69457216785935 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.5840124850865 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,24), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,24), colour="gray") +
  stat_function(fun=function(x){(0.55738362930159 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,24), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(-1.2075816360495 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-1.041156406882 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-1.179552391459 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.62417998266785 + 0.1502625498244*x + -0.007152664472829*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.22050751507995*x + -0.007351369134945*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.8029256146075 + 0.165957025626288*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.305072810695*x + -0.010255615653915*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.531896660533635 + 0.15096970501841*x + -0.0070283433490403*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.28232386593*x + -0.010855039139*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.1618843579945*x + -0.007749524154465*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.612258174813192 + 0.1624878116603*x + -0.00677339887530375*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.8765822027105 + 0.12783802081225*x + -0.0060584714911248*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.161764513989639*x + -0.00793768915954*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.60103746678025 + 0.144541498191466*x + -0.00707287306606*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.16419771266985*x + -0.00855055603015*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(0 + 0.18871553213095*x + -0.0079616198584991*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.6825435607292 + 0.137650032442615*x + -0.0069294617549395*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,25), colour="gray") +
  stat_function(fun=function(x){(-0.8489529958284 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.78119335871065 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.9037523690473 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.719734792237655 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(-0.82242632105255 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0.17502039186843*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,19), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,19), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,9), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,16), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,16), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,16), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(-0.64162229108745 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(-0.593067034201 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(-0.646311489052505 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(-0.5393776845379 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,12), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-0.954688543822425 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-0.8696473892294 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-0.79913355414625 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-1.223717799137 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-1.246553061853 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-1.073649414203 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-0.80168726255 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-0.97178867216685 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(-1.0034198361305 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(-1.207569756555 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(-0.978835861498 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(1.05995426871655 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(-0.555822251316005 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(-0.668371769321515 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(-0.6843305145485 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(-0.62748888184675 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(-0.6545688189035 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,6), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,11), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,8), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,4), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,10), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,27), colour="gray") +
  stat_function(fun=function(x){(-0.8106907867263 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,27), colour="gray") +
  stat_function(fun=function(x){(-0.756132353018845 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,27), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,18), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,7), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0.3369178605105*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.9030524820745 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.71381111444355 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.735794265635 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.5709567413762 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.6379321622963 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.77099010806735 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.846431544735 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.70597576568755 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(-0.7160050317985 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,5), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.6112685 + 0.362312}, size=0.5, xlim=c(0,3), colour="gray") +
  
  

  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.36728550 + 0*x + -0*x^2)*0.6112685 + 0.362312}, size=3, xlim=c(0,23), colour='black')
print(meanPlot)
#export at 1200x1000


###by resource mani

#mean change
meanResourceDrought <- mean%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
meanResource <- mean%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(meanResourceDrought)
  
# meanResourcePlot <- ggplot(data=barGraphStats(data=meanResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.5)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')

#dispersion change

dispersionResourceDrought <- dispersion%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
dispersionResource <- dispersion%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(dispersionResourceDrought)

# dispersionResourcePlot <- ggplot(data=barGraphStats(data=dispersionResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.08, 0.15)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.14, label='(b)', size=10, hjust='left')

#richness change

richnessResourceDrought <- richness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
richnessResource <- richness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(richnessResourceDrought)

# richnessResourcePlot <- ggplot(data=barGraphStats(data=richnessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.22, 0.16)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')

#evenness change
evennessResourceDrought <- evenness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
evennessResource <- evenness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(evennessResourceDrought)

# evennessResourcePlot <- ggplot(data=barGraphStats(data=evennessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.06)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left')
  
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))


#by resource at final year
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=meanResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(0, 0.5)) +
  xlab('')+
  annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')
dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=dispersionResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.08, 0.15)) +
  xlab('')+
  annotate('text', x=0.5, y=0.14, label='(b)', size=10, hjust='left')
richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=richnessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.22, 0.16)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')
evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=evennessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(0, 0.06)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))

# #20 year estimate
# meanResourcePlot20 <- ggplot(data=barGraphStats(data=meanResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.5)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')
# dispersionResourcePlot20 <- ggplot(data=barGraphStats(data=dispersionResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.08, 0.15)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.14, label='(b)', size=10, hjust='left')
# richnessResourcePlot20 <- ggplot(data=barGraphStats(data=richnessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.22, 0.16)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')
# evennessResourcePlot20 <- ggplot(data=barGraphStats(data=evennessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.06)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlot20, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlot20, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlot20, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlot20, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))




#by resource mani boxplot
meanResourceBoxFinal <- ggplot(data=meanResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(0, 1.0)) +
  xlab('')+
  annotate('text', x=0.5, y=0.98, label='(a)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))
dispersionResourceBoxFinal <- ggplot(data=dispersionResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.1), name='Dispersion Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.4, 0.4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.39, label='(b)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))
richnessResourceBoxFinal <- ggplot(data=richnessResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-1, 0.6, 0.3), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-1, 0.6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.58, label='(c)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))
evennessResourceBoxFinal <- ggplot(data=evennessResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-0.4, 0.6, 0.2), name='Evenness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.4, 0.6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.58, label='(d)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))


pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourceBoxFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourceBoxFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourceBoxFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourceBoxFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))



#other inset options

#inset - density plot of mean change (all datapoints)
mean10yr <- ggplot(data=mean, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
mean20yr <- ggplot(data=mean, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
meanFinalYear <- ggplot(data=mean, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')
dispersion10yr <- ggplot(data=dispersion, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
dispersion20yr <- ggplot(data=dispersion, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
dispersionFinalYear <- ggplot(data=dispersion, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')
richness10yr <- ggplot(data=richness, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
richness20yr <- ggplot(data=richness, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
richnessFinalYear <- ggplot(data=richness, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')
evenness10yr <- ggplot(data=evenness, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
evenness20yr <- ggplot(data=evenness, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
evennessFinalYear <- ggplot(data=evenness, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')


pushViewport(viewport(layout=grid.layout(3,4)))
print(mean10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(mean20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(meanFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(dispersion10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(dispersion20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(dispersionFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(richness10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(richness20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
print(richnessFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 3))
print(evenness10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(evenness20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
print(evennessFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 4))

      
      
      
###summary stats from bayesian output
chainsCommunitySummary <- chainsCommunity%>%
  select(U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4, E_int.1.1, E_int.2.1, E_int.3.1, E_int.4.1, E_slope.1.1, E_slope.2.1, E_slope.3.1, E_slope.4.1, E_quad.1.1, E_quad.2.1, E_quad.3.1, E_quad.4.1, E_int.1.2, E_int.2.2, E_int.3.2, E_int.4.2, E_slope.1.2, E_slope.2.2, E_slope.3.2, E_slope.4.2, E_quad.1.2, E_quad.2.2, E_quad.3.2, E_quad.4.2, D_int.1.1, D_int.2.1, D_int.3.1, D_int.4.1, D_slope.1.1, D_slope.2.1, D_slope.3.1, D_slope.4.1, D_quad.1.1, D_quad.2.1, D_quad.3.1, D_quad.4.1, D_int.1.2, D_int.2.2, D_int.3.2, D_int.4.2, D_slope.1.2, D_slope.2.2, D_slope.3.2, D_slope.4.2, D_quad.1.2, D_quad.2.2, D_quad.3.2, D_quad.4.2)%>%
  gather(key=parameter, value=value, U_int.2.1:D_quad.4.2)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(CI=sd*2)

write.csv(chainsCommunitySummary, 'bayesian_summary stats.csv')
chainsCommunitySummary2 <- read.csv('bayesian_summary stats.csv')

#mean plots
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

#dispersion plots
dispersionIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

dispersionSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

dispersionQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

#richness plots
richnessIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

richnessSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

richnessQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

#evenness plots
evennessIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

evennessSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

evennessQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

pushViewport(viewport(layout=grid.layout(4,3)))
print(meanIntPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanSlopePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanQuadPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(dispersionIntPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(dispersionSlopePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(dispersionQuadPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
print(richnessIntPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(richnessSlopePlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(richnessQuadPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 3))
print(evennessIntPlot, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(evennessSlopePlot, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(evennessQuadPlot, vp=viewport(layout.pos.row = 4, layout.pos.col = 3))
#export at 2400x2000


#overall responses
meanOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='mean change'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='dispersion'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='richness'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='evenness'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

pushViewport(viewport(layout=grid.layout(1,4)))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
#export at 1800x1200
