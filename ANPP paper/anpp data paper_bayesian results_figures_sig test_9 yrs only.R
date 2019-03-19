library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

#kim's laptop
setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

#kim's desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
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
#experiment information --------------------------------------------------------
expRaw <- read.csv('ExperimentInformation_Dec2016.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('ForBayesianAnalysisANPP_9yr_Dec2016.csv')

rawData2<- rawData%>%
  filter(plot_mani<6, anpp!='NA', anpp_PC!='NA')%>%
  summarise(mean=mean(anpp_PC), std=sd(anpp_PC)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani))

#for table of experiment summarizing various factors
expInfoSummary <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich), anpp=mean(anpp),
            MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  summarise(length_median=median(experiment_length), length_min=min(experiment_length), length_max=max(experiment_length),
            plot_mani_median=median(plot_mani), plot_mani_min=min(plot_mani), plot_mani_max=max(plot_mani),
            rrich_median=median(rrich), rrich_min=min(rrich), rrich_max=max(rrich),
            anpp_median=median(anpp), anpp_min=min(anpp), anpp_max=max(anpp),
            MAP_median=median(MAP), MAP_min=min(MAP), MAP_max=max(MAP),
            MAT_median=median(MAT), MAT_min=min(MAT), MAT_max=max(MAT))%>%
  gather(variable, estimate)

#treatment info
trtInfo <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich),
            anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(resource_mani=(nutrients+carbon+irrigation+drought), id=1:length(treatment))

################################################################################
################################################################################

# #only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# chains1 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_9yr\\anpp_9yr_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_9yr\\anpp_9yr_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_9yr\\anpp_9yr_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_9yr\\anpp_9yr_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)
# 
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsCommunity$mu.1))
# plot(density(chainsCommunity$mu.2))
# plot(density(chainsCommunity$mu.3))
# 
# 
# #get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
# chainsCommunity2 <- chainsCommunity%>%
#   select(lp__,
#          #plot_mani intercepts (first digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.1, U.2.1, U.3.1, U.4.1,
#          #plot_mani linear slopes (first digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.2, U.2.2, U.3.2, U.4.2,
#          #plot_mani quad slopes (first digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.3, U.2.3, U.3.3, U.4.3,
#          #ANPP intercept, linear, and quad slopes (first digit): 1=anpp
#          D.1.1, D.1.2, D.1.3,
#          #richness intercept, linear, and quad slopes (first digit): 2=richness
#          D.2.1, D.2.2, D.2.3,
#          #MAP intercept, linear, and quad slopes (first digit): 1=MAP
#          E.1.1, E.1.2, E.1.3,
#          #MAT intercept, linear, and quad slopes (first digit): 2=MAT
#          E.2.1, E.2.2, E.2.3,
#          #overall intercept, linear, and quad slopes
#          mu.1, mu.2, mu.3)%>%
#   gather(key=parameter, value=value, U.1.1:mu.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))
# 
# write.csv(chainsCommunity2, 'bayesian_output_summary_anpp_03132017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_anpp_03132017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,1779:2501]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 1779:2501])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,1779:2501]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 1779:2501])'] <- 'sd'
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
#   #spread by parameter
#   select(id, parameter, mean)%>%
#   spread(key=parameter, value=mean)
# 
# write.csv(chainsFinal, 'bayesian_output_mean sd_anpp_03132017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_anpp_03132017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  mutate(yr9=(intercept+linear*9+quadratic*9^2)*(0.516839)+(0.3145941))%>%
  mutate(yr_final=(intercept+linear*alt_length+quadratic*alt_length^2)*(0.516839)+(0.3145941))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_anpp.csv', row.names=F)



###main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  # coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('ANPP Change')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(1.04552526203964 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.97068517552708 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.753492279261244 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8635072428376 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.0519910325252 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.718400253719 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.74621284573844 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7125129034316 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.70603183861024 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.787994527102876 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.77525706737756 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.865525974827488 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.82218057543436 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.508966793564 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.920258884171248 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.90725968028448 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(1.83122928852 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-1.3108444386668 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.87024816559396 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(1.307429745896 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(2.29603408612 + -0.38599446725868*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(3.5798466744 + -0.41634272799172*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.57532283876664*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(2.7541002928 + -0.94965914968*x + 0.0767325030632*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(3.4174605988 + -1.00663111764*x + 0.0751723342644*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(1.17585863265 + -1.106867822816*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(3.5005715404 + -1.5041220980484*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(4.159653992 + -1.5731626062504*x + 0.1372080005772*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.83242814352 + -0.7232439532*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(2.5027853924 + -0.723041681408*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.582507247777292*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.91018114134444 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.39242254112816*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3856510837036*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.38722433705136*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.32295475821984*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8219430156908 + 0.38007864220196*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.4134937025524*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.40010487450432*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.64179898669162 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.678188113878292 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.79975029709052 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.75589940097122 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8614544688134 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.72526915424656 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.68344311231262 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.3896418729076 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.2613819503632 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.3043891728514 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.660812896177708 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6490209233048 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.039404904572 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1736705274372 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.91655541121584 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0744092295496 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.688456992832 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.7636535147704 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.609001340644 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.912519415677076 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.4316886834268 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.958154322121608 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1538248955024 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.655104014104 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0785620284356 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(1.876129916132 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.686855937084366 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.747757967981314 + 0.6542416914468*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.26503937502 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.0822450312188 + 0*x + -0.076922960757488*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.44359216434964*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.08658571682 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0554191765064 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines (plot mani 1-4 staggered by intercept so lines don't overlap)
  #mani1
  stat_function(fun=function(x){(-0.40832900 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.40832900+0.27062850) + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){((-0.40832900+0.49680000) + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.40832900 + 0*x + 0*x^2)*(0.516839)+(0.3145941)}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.40832900+2.59186500) + (0.120001500-0.66454400)*x + (-0.016923000+0.07270935)*x^2)*(0.516839)+(0.3145941)}, size=3, xlim=c(0,8), colour='#EC1804')

print(meanPlot) #export at 1200x1000


##summary stats from bayesian output --------------------------------------------------------
#gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(#plot_mani intercepts (first digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.1, U.2.1, U.3.1, U.4.1,
#          #plot_mani linear slopes (first digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.2, U.2.2, U.3.2, U.4.2,
#          #plot_mani quad slopes (first digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.3, U.2.3, U.3.3, U.4.3,
#          #ANPP intercept, linear, and quad slopes (first digit): 1=anpp
#          D.1.1, D.1.2, D.1.3,
#          #richness intercept, linear, and quad slopes (first digit): 2=richness
#          D.2.1, D.2.2, D.2.3,
#          #MAP intercept, linear, and quad slopes (first digit): 1=MAP
#          E.1.1, E.1.2, E.1.3,
#          #MAT intercept, linear, and quad slopes (first digit): 2=MAT
#          E.2.1, E.2.2, E.2.3,
#          #overall intercept, linear, and quad slopes
#          mu.1, mu.2, mu.3)%>%
#   gather(key=parameter, value=value, U.1.1:mu.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(CI=sd*2)%>%
#   separate(parameter, c('level', 'predictor', 'parameter'))%>%
#   mutate(parameter=ifelse(level=='mu', predictor, parameter), predictor=ifelse(level=='mu', 'overall', predictor))%>%
#   #rename parts to be more clear
#   mutate(parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          predictor=ifelse(level=='D'&predictor==1, 'ANPP', ifelse(level=='D'&predictor==2, 'rrich', ifelse(level=='E'&predictor==1, 'MAP', ifelse(level=='E'&predictor==2, 'MAT', ifelse(level=='U'&predictor==1, 'plot mani 2', ifelse(level=='U'&predictor==2, 'plot mani 3', ifelse(level=='U'&predictor==3, 'plot mani 4', ifelse(level=='U'&predictor==4, 'plot mani 5', 'overall')))))))))%>%
#   select(level, parameter, predictor, median, sd, CI)
# 
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_anpp_03132017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_anpp_03132017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)

###overall responses from bayesian output --------------------------------------------------------
ggplot(data=subset(chainsCommunityOverall, predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('ANPP Change')
#export at 800x800

#other predictor plots --------------------------------------------------------
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.15, 3.5) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 0.5) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.15, 0.1), limits=c(-0.13,0.16)) +
  coord_flip()

#plot all together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(3,1))) 
print(meanIntPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanSlopePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(meanQuadPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
#export at 800x1200



###by resource mani - raw data--------------------------------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

rawTrt <- rawData%>%
  filter(anpp_PC!='NA')%>%
  filter(treatment_year==9|treatment_year==experiment_length)%>%
  select(site_code, project_name, community_type, treatment, plot_mani, rrich, anpp, MAT, MAP, experiment_length, treatment_year, anpp_PC)%>%
  left_join(trtDetail)%>%
  mutate(resource_mani=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
                              ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
                                     ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
                                            ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
                                                   ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'nuts:CO2',
                                                          ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'nuts:dro',
                                                                 ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'nuts:irr',
                                                                        ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'CO2:dro',
                                                                               ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'CO2:irr',
                                                                                      ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'nuts:CO2:dro', 
                                                                                             ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr','other'))))))))))))%>%
  mutate(resource_mani_combo=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
                                    ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
                                           ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
                                                  ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
                                                         ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'multiple',
                                                                ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'multiple',
                                                                       ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'multiple',
                                                                              ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'multiple',
                                                                                     ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'multiple',
                                                                                            ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'multiple', 
                                                                                                   ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'multiple','other'))))))))))))

#plot raw data by resource manipulated at final year of each experiment (varies by experiment) ---------------------------
ggplot(data=barGraphStats(data=rawTrt, variable='anpp_PC', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.10), name='ANPP Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.4, 0.7), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.7, label='(a)', size=12, hjust='left') +
  annotate('text', x=1, y=0.12, label='a', size=10) +
  annotate('text', x=2, y=0.51, label='b*', size=10) +
  annotate('text', x=3, y=0.18, label='c*', size=10) +
  annotate('text', x=4, y=0.35, label='abc*', size=10) +
  annotate('text', x=5, y=-0.36, label='d*', size=10) +
  annotate('text', x=6, y=0.68, label='abc', size=10)
#export at 1200 x 1000


###by magnitude of resource manipulated---------------------------------
Nplot <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=anpp_PC)) +
  geom_point(size=5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-3, 2.5, 0.50), name='ANPP Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.9, y=2.5, label='(a)', size=12, hjust='left')
H2Oplot <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=anpp_PC)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(-3, 3, 0.50), name='') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-51, y=1.5, label='(b)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(1,2)))
print(Nplot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(H2Oplot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2000 x 1000