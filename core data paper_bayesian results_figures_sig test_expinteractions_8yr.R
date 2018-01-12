library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)
library(nlme)

#kim's laptop
setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

#kim's desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
             plot.title = element_text(size=40, vjust=2),
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
expRaw <- read.csv('ExperimentInformation_Nov2017.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

trtInfo <- read.csv('ExperimentInformation_Nov2017.csv')%>%
  select(-X)

rawData <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\8yr or less\\ForAnalysis_allAnalysisShort.csv')

rawData2<- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_rich=mean(S_PC), std_rich=sd(S_PC)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length))

#for table of experiment summarizing various factors
expInfoSummary <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich), anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  summarise(length_mean=mean(experiment_length), length_min=min(experiment_length), length_max=max(experiment_length),
            plot_mani_median=mean(plot_mani), plot_mani_min=min(plot_mani), plot_mani_max=max(plot_mani),
            rrich_mean=mean(rrich), rrich_min=min(rrich), rrich_max=max(rrich),
            anpp_mean=mean(anpp), anpp_min=min(anpp), anpp_max=max(anpp),
            MAP_mean=mean(MAP), MAP_min=min(MAP), MAP_max=max(MAP),
            MAT_mean=mean(MAT), MAT_min=min(MAT), MAT_max=max(MAT))%>%
  gather(variable, estimate)

#treatment info
trtInfo <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, trt_type)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich),
            anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(resource_mani=(nutrients+carbon+irrigation+drought), id=1:length(treatment))

studyInfo <- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  select(site_code, community_type, project_name, treatment)%>%
  unique()


################################################################################
################################################################################

# #only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\8yr or less\\simple_mod_effs_coding_interactions_short_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\8yr or less\\simple_mod_effs_coding_interactions_short_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\8yr or less\\simple_mod_effs_coding_interactions_short_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\8yr or less\\simple_mod_effs_coding_interactions_short_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)
# 
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsCommunity$D.1.1.1))
# plot(density(chainsCommunity$D.1.1.2))
# plot(density(chainsCommunity$D.1.1.3))
# 
# 
# #get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
# #mean change are the 1's, richness are the 2's
# chainsCommunity2 <- chainsCommunity%>%
#   select(lp__,
#          #trt_type intercepts: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.1, U.2.1.1, U.1.2.1, U.2.2.1, U.1.3.1, U.2.3.1, U.1.4.1, U.2.4.1, U.1.5.1, U.2.5.1,
#          U.1.6.1, U.2.6.1, U.1.7.1, U.2.7.1, U.1.8.1, U.2.8.1, U.1.9.1, U.2.9.1, U.1.10.1, U.2.10.1,
#          U.1.11.1, U.2.11.1, U.1.12.1, U.2.12.1, U.1.13.1, U.2.13.1, U.1.14.1, U.2.14.1, U.1.15.1, U.2.15.1,
#          U.1.16.1, U.2.16.1, U.1.17.1, U.2.17.1, U.1.18.1, U.2.18.1, U.1.19.1, U.2.19.1, U.1.20.1, U.2.20.1,
#          U.1.21.1, U.2.21.1, U.1.22.1, U.2.22.1, U.1.23.1, U.2.23.1, U.1.24.1, U.2.24.1, U.1.25.1, U.2.25.1,
#          U.1.26.1, U.2.26.1, U.1.27.1, U.2.27.1, U.1.28.1, U.2.28.1, U.1.29.1, U.2.29.1, U.1.30.1, U.2.30.1,
#          U.1.31.1, U.2.31.1, U.1.32.1, U.2.32.1, U.1.33.1, U.2.33.1, U.1.34.1, U.2.34.1, U.1.35.1, U.2.35.1,
#          U.1.36.1, U.2.36.1, U.1.37.1, U.2.37.1, U.1.38.1, U.2.38.1, U.1.39.1, U.2.39.1, U.1.40.1, U.2.40.1,
#          U.1.41.1, U.2.41.1, U.1.42.1, U.2.42.1, U.1.43.1, U.2.43.1, U.1.44.1, U.2.44.1, U.1.45.1, U.2.45.1,
#          U.1.46.1, U.2.46.1, U.1.47.1, U.2.47.1, U.1.48.1, U.2.48.1, U.1.49.1, U.2.49.1, U.1.50.1, U.2.50.1,
#          U.1.51.1, U.2.51.1, U.1.52.1, U.2.52.1, U.1.53.1, U.2.53.1, U.1.54.1, U.2.54.1,
#          #trt_type linear slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.2, U.2.1.2, U.1.2.2, U.2.2.2, U.1.3.2, U.2.3.2, U.1.4.2, U.2.4.2, U.1.5.2, U.2.5.2,
#          U.1.6.2, U.2.6.2, U.1.7.2, U.2.7.2, U.1.8.2, U.2.8.2, U.1.9.2, U.2.9.2, U.1.10.2, U.2.10.2,
#          U.1.11.2, U.2.11.2, U.1.12.2, U.2.12.2, U.1.13.2, U.2.13.2, U.1.14.2, U.2.14.2, U.1.15.2, U.2.15.2,
#          U.1.16.2, U.2.16.2, U.1.17.2, U.2.17.2, U.1.18.2, U.2.18.2, U.1.19.2, U.2.19.2, U.1.20.2, U.2.20.2,
#          U.1.21.2, U.2.21.2, U.1.22.2, U.2.22.2, U.1.23.2, U.2.23.2, U.1.24.2, U.2.24.2, U.1.25.2, U.2.25.2,
#          U.1.26.2, U.2.26.2, U.1.27.2, U.2.27.2, U.1.28.2, U.2.28.2, U.1.29.2, U.2.29.2, U.1.30.2, U.2.30.2,
#          U.1.31.2, U.2.31.2, U.1.32.2, U.2.32.2, U.1.33.2, U.2.33.2, U.1.34.2, U.2.34.2, U.1.35.2, U.2.35.2,
#          U.1.36.2, U.2.36.2, U.1.37.2, U.2.37.2, U.1.38.2, U.2.38.2, U.1.39.2, U.2.39.2, U.1.40.2, U.2.40.2,
#          U.1.41.2, U.2.41.2, U.1.42.2, U.2.42.2, U.1.43.2, U.2.43.2, U.1.44.2, U.2.44.2, U.1.45.2, U.2.45.2,
#          U.1.46.2, U.2.46.2, U.1.47.2, U.2.47.2, U.1.48.2, U.2.48.2, U.1.49.2, U.2.49.2, U.1.50.2, U.2.50.2,
#          U.1.51.2, U.2.51.2, U.1.52.2, U.2.52.2, U.1.53.2, U.2.53.2, U.1.54.2, U.2.54.2,
#          #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.3, U.2.1.3, U.1.2.3, U.2.2.3, U.1.3.3, U.2.3.3, U.1.4.3, U.2.4.3, U.1.5.3, U.2.5.3,
#          U.1.6.3, U.2.6.3, U.1.7.3, U.2.7.3, U.1.8.3, U.2.8.3, U.1.9.3, U.2.9.3, U.1.10.3, U.2.10.3,
#          U.1.11.3, U.2.11.3, U.1.12.3, U.2.12.3, U.1.13.3, U.2.13.3, U.1.14.3, U.2.14.3, U.1.15.3, U.2.15.3,
#          U.1.16.3, U.2.16.3, U.1.17.3, U.2.17.3, U.1.18.3, U.2.18.3, U.1.19.3, U.2.19.3, U.1.20.3, U.2.20.3,
#          U.1.21.3, U.2.21.3, U.1.22.3, U.2.22.3, U.1.23.3, U.2.23.3, U.1.24.3, U.2.24.3, U.1.25.3, U.2.25.3,
#          U.1.26.3, U.2.26.3, U.1.27.3, U.2.27.3, U.1.28.3, U.2.28.3, U.1.29.3, U.2.29.3, U.1.30.3, U.2.30.3,
#          U.1.31.3, U.2.31.3, U.1.32.3, U.2.32.3, U.1.33.3, U.2.33.3, U.1.34.3, U.2.34.3, U.1.35.3, U.2.35.3,
#          U.1.36.3, U.2.36.3, U.1.37.3, U.2.37.3, U.1.38.3, U.2.38.3, U.1.39.3, U.2.39.3, U.1.40.3, U.2.40.3,
#          U.1.41.3, U.2.41.3, U.1.42.3, U.2.42.3, U.1.43.3, U.2.43.3, U.1.44.3, U.2.44.3, U.1.45.3, U.2.45.3,
#          U.1.46.3, U.2.46.3, U.1.47.3, U.2.47.3, U.1.48.3, U.2.48.3, U.1.49.3, U.2.49.3, U.1.50.3, U.2.50.3,
#          U.1.51.3, U.2.51.3, U.1.52.3, U.2.52.3, U.1.53.3, U.2.53.3, U.1.54.3, U.2.54.3,
#          #ANPP intercept, linear, and quad slopes (center digit): 2=anpp
#          D.1.2.1, D.2.2.1,
#          D.1.2.2, D.2.2.2,
#          D.1.2.3, D.2.2.3,
#          #richness intercept, linear, and quad slopes (center digit): 3=gamma diversity
#          D.1.3.1, D.2.3.1,
#          D.1.3.2, D.2.3.2,
#          D.1.3.3, D.2.3.3,
#          #overall intercept, linear, and quad slopes (center digit): 1=overall
#          D.1.1.1, D.2.1.1,
#          D.1.1.2, D.2.1.2,
#          D.1.1.3, D.2.1.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:D.2.1.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))
# 
# write.csv(chainsCommunity2, 'bayesian_output_summary_expinteractions_short_01122018.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteractions_short_01122018.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,3158:5377]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 3158:5377])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,3158:5377]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 3158:5377])'] <- 'sd'
# 
# chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
#   #split names into parts
#   separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
#   select(-B)%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', 'richness'),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          id=as.integer(id))%>%
#   #if 95% confidence interval overlaps 0, then set mean to 0
#   mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
#   #spread by variable
#   select(variable, id, parameter, mean)%>%
#   spread(key=parameter, value=mean)
# 
# write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_short_01122018.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_short_01122018.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=8, 7, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1404912)+(0.2887009), (intercept+linear*7+quadratic*7^2)*(0.2275677)+(-0.03003684)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1404912)+(0.2887009),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2275677)+(-0.03003684)))%>%
  mutate(color=ifelse(anpp<281&variable=='mean', '#1104DC44', ifelse(anpp<561&anpp>280&variable=='mean', '#4403AE55', ifelse(anpp<841&anpp>560&variable=='mean', '#77038166', ifelse(anpp>840&variable=='mean', '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,',
                       '*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,'),
         curve5='), colour=',
         curve6=') +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_short_01122018.csv', row.names=F)

#summary lines
chainsEquationsSummary <- chainsEquations%>%
  group_by(variable, trt_overall)%>%
  summarise(intercept_mean=mean(intercept), intercept_sd=sd(intercept), linear_mean=mean(linear), linear_sd=sd(linear), quadratic_mean=mean(quadratic), quadratic_sd=sd(quadratic))%>%
  ungroup()%>%
  mutate(intercept_high=intercept_mean+1.96*intercept_sd, intercept_low=intercept_mean-1.96*intercept_sd, linear_high=linear_mean+1.96*linear_sd, linear_low=linear_mean-1.96*linear_sd, quadratic_high=quadratic_mean+1.96*quadratic_sd, quadratic_low=quadratic_mean-1.96*quadratic_sd)%>%
  mutate(intercept_high2=ifelse(intercept_high>0, 1, 0), intercept_low2=ifelse(intercept_low<0, 1, 0), linear_high2=ifelse(linear_high>0, 1, 0), linear_low2=ifelse(linear_low<0, 1, 0), quadratic_high2=ifelse(quadratic_high>0, 1, 0), quadratic_low2=ifelse(quadratic_low<0, 1, 0))%>%
  mutate(intercept_cross=intercept_high2+intercept_low2, linear_cross=linear_high2+linear_low2, quadratic_cross=quadratic_high2+quadratic_low2)%>%
  mutate(intercept_final=ifelse(intercept_cross==2, 0, intercept), linear_final=ifelse(linear_cross==2, 0, linear), quadratic_final=ifelse(quadratic_cross==2, 0, quadratic))

###main figure (Figure 1)
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1), labels=seq(1,8,1)) +
  ylim(-10,10) +
xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  annotate('text', x=0, y=1, label='(b)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.54625653153552 + 0.43972066798595*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7859746237845 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0.52307253683465 + 0.35139291315061*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0.77335640939415 + 0.34583986589495*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.4181501867242*x + -0.049934832944815*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8167570661265 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.59388272241365 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.01294972944 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.76575863656805 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.840697515615 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.937475036838 + 0.55338008641*x + -0.0670897978936*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.4361330661 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7041445445815 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.49703410545 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.33471045765 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.22179114695 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.490190594476*x + -0.052986409480152*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.848109171874 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.99992687575 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8581998848455 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8259338522875 + 0.395891730620975*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0087025232175 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.5950985945185*x + -0.0692312085281*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.842384938965 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7353055091505 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.733460846116 + 0.44619816785715*x + -0.085449115686935*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9070845916635 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.67936734494192 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.9133193168195 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.175097351288 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.0521322819 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.3328541 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.217312235395 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.1779233243 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.801861931955 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.57822852947435 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.980141192453 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.80092888264335 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7210495590047 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.818109113681 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.8091666481685 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.62277671403175 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.836675988323593 + 0.5098724240965*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.002122608465 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6013024734035 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.282086897300195*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4512384490549*x + -0.04390282565553*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.373318489423751*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.68332573416165 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7121351943654 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.8238312050615 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.01807942092 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.0687765398 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.1808183036 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.400849670253645*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.9337251183825 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.54381163863185 + 0.38056163784405*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.62095718428495 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0.85168894414855 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.87773240514*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.560215246119*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.62189554156325 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.589154177562 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.850057178616 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.791221264736 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.16971832875 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6169498731345 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.175669139 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.28224997685 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.30898565915 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.076284905752 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.804528490475 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7130035773346 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.36815918085 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6199611242501 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.728007044677 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6402577418467 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.61029349029995 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.336268066589495*x + -0.06447664235354*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.707490387025 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7890085091472 + 0.52643002353625*x + -0.0787453280632*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.30714778695 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.14221421195 + 0.420594001572*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.00088704342 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.99519312494 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7264198903645 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.77813288848 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.4799449734243*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.06095321405 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.529255745940925 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.471034111326035*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.46008884015825*x + -0.084269049484141*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.715156444875 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8045124149145 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.64403337018625 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.61410942021105 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7873501024755 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.40665457035825*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6104036246015 + 0.36119715621085*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.807772481098 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.878529289576 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.30272229235 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.086019159025 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.63919718106805 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.821876404598 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.98322530175 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.58765116811248 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.781419799514 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.757582848256 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.8823349573795 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.577311488229585 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7015496751561 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.75587219106795 + 0*x + 0.049706389068421*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7143462939895 + 0*x + 0.0518018509909*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.952003337981345 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8572499295015 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + -0.10511772538602*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4644904297168*x + -0.0866228062048785*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7902363251285 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.52358182279195*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.40559345412625*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.36914446126945*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.60621666560245 + 0.42706301221395*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0.5247156836514 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0.495973394084 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.742008130292*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0.56584655672625 + 0.67053477956*x + -0.074393053337045*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.42787614715375*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.53381092702845*x + -0.05928521760251*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.21086100265 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.8216386384665 + 0.4831739154995*x + -0.06598802795255*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.67974312656755 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.880520813021 + 0.526919309318*x + -0.060120187448775*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.3684279703 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.46326359245 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.3313588766 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.635671747433 + 0.4131379532197*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.54421001994825 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6055534886457 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7975137264705 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.749388387427 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.02904190928 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.09496870603 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.221488980255 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.885895151335 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.82626783354727 + 0.4657021965833*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.59202228292815 + 0.28247474340477*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.52300955639765 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.51140721309015 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.4980076780983*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.5425147684863*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.40597228586125*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.5617447969848*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.49440596935534 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.667753626820075 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7014897078679 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-1.0031347394 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.54579104326*x + -0.05172593560662*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.51807125321795*x + -0.0585348791316*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.501515529295*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9124160350625 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.33030176635865*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.617633943805*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.52053809768*x + -0.049120395110599*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6387075182127 + 0.5633910117565*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.0919549274075 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.298166271 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.1294238772 + 0.3329089965046*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.1742604071 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0.65766107536615 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + -1.03876744409325*x + 0.4696892133972*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + -1.2433620178905*x + 0.49695869022075*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0936402802329995*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0.59056053856655 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.63418674278845 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.895628177346 + 0.45957135197885*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.2500562438 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.06575692426 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.96208572677 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.05940799565 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(-0.7258152971405 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(-0.6523586991233 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.48528248447245*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0.70805279581185 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.5905819386 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.97145436818 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.87868476898 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.495922403243*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.460825782045775*x + -0.087952840992015*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(1.26841639737 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.47052978445445*x + -0.07787852275505*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.64662620458285 + 0.4041046562996*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.46355033719928*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.51074635545 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.70290119855 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.47267210815 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.63752942645 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.71409253655 + 0.49604033242035*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.5598986463 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.61212185038505 + 0.3369717706043*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.7501841341*x + -0.069046909055*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0.55273997783325 + 0.83513568695*x + -0.08456857361665*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8472963792448 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.763879170601 + 0.4045085877174*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.72207292012 + 0.6456994011738*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6668786317478 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 1.13385133415*x + -0.12163120252*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.70998512512*x + -0.05259614371692*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.5859016452315*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.7593047031*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(-0.73304483660195 + 0.685612463305*x + -0.06677491747825*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.33472008960925*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.94474384635*x + -0.05486814193374*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0.8850705357255 + 1.073325381*x + -0.07957527412115*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.574426401613065 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.54963658074382 + 0.522904212115*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.88055266305*x + -0.05338500576934*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.069159320715 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0242137917085 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7991418820325 + 0.56073625246535*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.30928325079545*x + -0.05005948010845*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.04998602371339*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.04548876583447*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + -0.4633618447722*x + 0.0885066752673*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4816503595371*x + -0.10465244701481*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.3501463356054*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.668658636945*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.3829890872357*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.5801472448255*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.392310883439*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.916686199584555 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7092773905469 + 0.448884611041*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.067615897445 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.57739051926*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.725819034323775 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.01326174812 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.40935162211365*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.573472797136*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.69015387204*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.675696120843*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.894243373065 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.77545985730585 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0.809522405302 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.5973351463891 + 0.61436102514011*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.5734242751334*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1404912)+(0.2887009)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.4388275 + 0.2006625*x + -0.00000000*x^2)*(0.1404912)+(0.2887009)}, size=5, xlim=c(0,7), colour='black')

# print(meanPlot) #export at 1200x1000



#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,1.0))  +
  scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1), labels=seq(1,8,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('') +
  ylab('Richness Difference') +
  annotate('text', x=0, y=1.0, label='(a)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9129803024635 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.770915047794 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.92259858786275 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7990047935285 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.62974679473185 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.91344065628045 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.621278715711*x + -0.09461518709805*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8654698208065 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8821976630025 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.72305140611297 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.107752093499 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.137344005335 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9126683847489 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.80562323595 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.266051405148 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.87501526322295 + 1.458245420353*x + -0.47935304281245*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.82801783331025*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.0464355364694 + -0.70573570126635*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.76797814874215*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.70929396729145 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.077026922435 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.68920920264983 + 0.50033495798885*x + -0.069516006725533*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8530853448702 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.70244589553395 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.897178108287*x + -0.137969209352395*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.66560190409746*x + -0.109435801537195*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.827153163193865 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.356620452335 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.95524583366785 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.91664316329665 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.5720700237 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.75432705981369 + 0.71288877919175*x + -0.08680975087179*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 1.2043412154685*x + -0.3377508022588*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.72246875639725 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6011625158514 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.57629556767393*x + -0.090518569183082*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 1.05144525915*x + -0.14282964788*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.97362123*x + -0.12865711528*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7357771011918 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 1.0070334134859*x + -0.139825491162795*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.68629560073331 + -0.55518900499*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.66520455594055 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.674433120334845 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9792696433931 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.9609033084521 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.456119926350425*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.499787942138835*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.176681220305 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.9845689868864 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.280997555846 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.07399506348265 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.84201733366544 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.967163145024715 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.764273295704795*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.793044839593*x + -0.10309066142752*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.71917680857165*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6353438759872*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.792697700149585 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.53698363397586*x + -0.109950992752*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6458182421215*x + -0.124697069125*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.8893138278*x + -0.1232773962108*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 1.00412043075*x + -0.143488591865*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.87648390681*x + -0.1135907453858*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.861467639*x + -0.1104062592355*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 1.262964845*x + -0.193630562735*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.824996196964 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.86827882305 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.06522540220145 + 1.0483522953205*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.648085591891975*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.569495752884345*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.65851892584215 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.75703919779465 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.07201580356102*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.71465853086358 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.22475888687115*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5905764906218*x + -0.099526861892617*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 2.29360814285*x + -0.36697315123945*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 1.0733232079469*x + -0.156866406829945*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.676084336153675 + -0.5051517512689*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -1.28950608555*x + 0.1207941596775*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -1.5801635029*x + 0.158531512285*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7618688569655 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.179346107425 + -0.9117448777825*x + 0.07374191607775*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -1.08562381097*x + 0.10042976312095*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.944944771105*x + 0.08980558734489*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.98799656115*x + 0.10325630808735*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9159990785345 + -0.889263508595*x + 0.10209639168385*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.9742763564*x + 0.113233154461*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.208279528*x + 0.099256451645*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -1.4095813842*x + 0.150047255905*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.43801733024081*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -1.20040067655*x + 0.118153345353*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -1.2724870503*x + 0.1235505975885*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.683298831768*x + 0.0599615197600195*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.623100373319075 + -0.86972998253*x + 0.08716484238745*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0630152319914965*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0628527110827695*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.088675404379315*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.447623667246*x + -0.075283360354495*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.495322784116006*x + -0.06903505517808*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.55563231075555*x + -0.08711992337105*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 1.371381250475*x + -0.4420906290409*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.082522476555288*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(1.04616651200485 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.84594470397636 + 1.38465941965*x + -0.17705757597*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 1.10286973125*x + -0.177351106665*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.9419566667095*x + -0.146965022575*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.95382262082435 + 1.25687701655*x + -0.19299495233*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 1.163240556575*x + -0.18282515435*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.497102512116015*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.04626638334 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2275677)+(-0.03003684)}, size=2, xlim=c(0,2), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.2275677 + -0.03003684}, size=5, xlim=c(0,7), colour='black')

# print(richnessPlot) #export at 1200x1000

#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,1)))
print(richnessPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(meanPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1200 x 2400



# #summary stats from bayesian output --------------------------------------------------------
# #gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(
#          #trt_type intercepts: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.1, U.2.1.1, U.1.2.1, U.2.2.1, U.1.3.1, U.2.3.1, U.1.4.1, U.2.4.1, U.1.5.1, U.2.5.1,
#          U.1.6.1, U.2.6.1, U.1.7.1, U.2.7.1, U.1.8.1, U.2.8.1, U.1.9.1, U.2.9.1, U.1.10.1, U.2.10.1,
#          U.1.11.1, U.2.11.1, U.1.12.1, U.2.12.1, U.1.13.1, U.2.13.1, U.1.14.1, U.2.14.1, U.1.15.1, U.2.15.1,
#          U.1.16.1, U.2.16.1, U.1.17.1, U.2.17.1, U.1.18.1, U.2.18.1, U.1.19.1, U.2.19.1, U.1.20.1, U.2.20.1,
#          U.1.21.1, U.2.21.1, U.1.22.1, U.2.22.1, U.1.23.1, U.2.23.1, U.1.24.1, U.2.24.1, U.1.25.1, U.2.25.1,
#          U.1.26.1, U.2.26.1, U.1.27.1, U.2.27.1, U.1.28.1, U.2.28.1, U.1.29.1, U.2.29.1, U.1.30.1, U.2.30.1,
#          U.1.31.1, U.2.31.1, U.1.32.1, U.2.32.1, U.1.33.1, U.2.33.1, U.1.34.1, U.2.34.1, U.1.35.1, U.2.35.1,
#          U.1.36.1, U.2.36.1, U.1.37.1, U.2.37.1, U.1.38.1, U.2.38.1, U.1.39.1, U.2.39.1, U.1.40.1, U.2.40.1,
#          U.1.41.1, U.2.41.1, U.1.42.1, U.2.42.1, U.1.43.1, U.2.43.1, U.1.44.1, U.2.44.1, U.1.45.1, U.2.45.1,
#          U.1.46.1, U.2.46.1, U.1.47.1, U.2.47.1, U.1.48.1, U.2.48.1, U.1.49.1, U.2.49.1, U.1.50.1, U.2.50.1,
#          U.1.51.1, U.2.51.1, U.1.52.1, U.2.52.1, U.1.53.1, U.2.53.1, U.1.54.1, U.2.54.1,
#          #trt_type linear slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.2, U.2.1.2, U.1.2.2, U.2.2.2, U.1.3.2, U.2.3.2, U.1.4.2, U.2.4.2, U.1.5.2, U.2.5.2,
#          U.1.6.2, U.2.6.2, U.1.7.2, U.2.7.2, U.1.8.2, U.2.8.2, U.1.9.2, U.2.9.2, U.1.10.2, U.2.10.2,
#          U.1.11.2, U.2.11.2, U.1.12.2, U.2.12.2, U.1.13.2, U.2.13.2, U.1.14.2, U.2.14.2, U.1.15.2, U.2.15.2,
#          U.1.16.2, U.2.16.2, U.1.17.2, U.2.17.2, U.1.18.2, U.2.18.2, U.1.19.2, U.2.19.2, U.1.20.2, U.2.20.2,
#          U.1.21.2, U.2.21.2, U.1.22.2, U.2.22.2, U.1.23.2, U.2.23.2, U.1.24.2, U.2.24.2, U.1.25.2, U.2.25.2,
#          U.1.26.2, U.2.26.2, U.1.27.2, U.2.27.2, U.1.28.2, U.2.28.2, U.1.29.2, U.2.29.2, U.1.30.2, U.2.30.2,
#          U.1.31.2, U.2.31.2, U.1.32.2, U.2.32.2, U.1.33.2, U.2.33.2, U.1.34.2, U.2.34.2, U.1.35.2, U.2.35.2,
#          U.1.36.2, U.2.36.2, U.1.37.2, U.2.37.2, U.1.38.2, U.2.38.2, U.1.39.2, U.2.39.2, U.1.40.2, U.2.40.2,
#          U.1.41.2, U.2.41.2, U.1.42.2, U.2.42.2, U.1.43.2, U.2.43.2, U.1.44.2, U.2.44.2, U.1.45.2, U.2.45.2,
#          U.1.46.2, U.2.46.2, U.1.47.2, U.2.47.2, U.1.48.2, U.2.48.2, U.1.49.2, U.2.49.2, U.1.50.2, U.2.50.2,
#          U.1.51.2, U.2.51.2, U.1.52.2, U.2.52.2, U.1.53.2, U.2.53.2, U.1.54.2, U.2.54.2,
#          #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.3, U.2.1.3, U.1.2.3, U.2.2.3, U.1.3.3, U.2.3.3, U.1.4.3, U.2.4.3, U.1.5.3, U.2.5.3,
#          U.1.6.3, U.2.6.3, U.1.7.3, U.2.7.3, U.1.8.3, U.2.8.3, U.1.9.3, U.2.9.3, U.1.10.3, U.2.10.3,
#          U.1.11.3, U.2.11.3, U.1.12.3, U.2.12.3, U.1.13.3, U.2.13.3, U.1.14.3, U.2.14.3, U.1.15.3, U.2.15.3,
#          U.1.16.3, U.2.16.3, U.1.17.3, U.2.17.3, U.1.18.3, U.2.18.3, U.1.19.3, U.2.19.3, U.1.20.3, U.2.20.3,
#          U.1.21.3, U.2.21.3, U.1.22.3, U.2.22.3, U.1.23.3, U.2.23.3, U.1.24.3, U.2.24.3, U.1.25.3, U.2.25.3,
#          U.1.26.3, U.2.26.3, U.1.27.3, U.2.27.3, U.1.28.3, U.2.28.3, U.1.29.3, U.2.29.3, U.1.30.3, U.2.30.3,
#          U.1.31.3, U.2.31.3, U.1.32.3, U.2.32.3, U.1.33.3, U.2.33.3, U.1.34.3, U.2.34.3, U.1.35.3, U.2.35.3,
#          U.1.36.3, U.2.36.3, U.1.37.3, U.2.37.3, U.1.38.3, U.2.38.3, U.1.39.3, U.2.39.3, U.1.40.3, U.2.40.3,
#          U.1.41.3, U.2.41.3, U.1.42.3, U.2.42.3, U.1.43.3, U.2.43.3, U.1.44.3, U.2.44.3, U.1.45.3, U.2.45.3,
#          U.1.46.3, U.2.46.3, U.1.47.3, U.2.47.3, U.1.48.3, U.2.48.3, U.1.49.3, U.2.49.3, U.1.50.3, U.2.50.3,
#          U.1.51.3, U.2.51.3, U.1.52.3, U.2.52.3, U.1.53.3, U.2.53.3, U.1.54.3, U.2.54.3,
#          #ANPP intercept, linear, and quad slopes (center digit): 2=anpp
#          D.1.2.1, D.2.2.1,
#          D.1.2.2, D.2.2.2,
#          D.1.2.3, D.2.2.3,
#          #richness intercept, linear, and quad slopes (center digit): 3=gamma diversity
#          D.1.3.1, D.2.3.1,
#          D.1.3.2, D.2.3.2,
#          D.1.3.3, D.2.3.3,
#          #overall intercept, linear, and quad slopes (center digit): 1=overall
#          D.1.1.1, D.2.1.1,
#          D.1.1.2, D.2.1.2,
#          D.1.1.3, D.2.1.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:D.2.1.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(CI=sd*2)%>%
#   separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', 'richness'),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          predictor2=ifelse(level=='D'&predictor==2, 'ANPP', ifelse(level=='D'&predictor==3, 'rrich', ifelse(level=='D'&predictor==1, 'overall', 'trt_type'))))%>%
#   select(level, parameter, variable, predictor, predictor2, median, sd, CI)
# 
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_short_01122018.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_short_01122018.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  mutate(type=paste(predictor2, parameter, sep='_'))



###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('rrich_quadratic', 'ANPP_quadratic', 'overall_quadratic', 'rrich_linear', 'ANPP_linear', 'overall_linear', 'rrich_intercept', 'ANPP_intercept', 'overall_intercept'),
                   labels=c('Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=3.5), linetype='dashed') +
  geom_vline(aes(xintercept=6.5), linetype='dashed') +
  coord_flip() +
  ggtitle('Community Difference') +
  annotate('text', x=9.2, y=-0.8, label='(a)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='richness' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('rrich_quadratic', 'ANPP_quadratic', 'overall_quadratic', 'rrich_linear', 'ANPP_linear', 'overall_linear', 'rrich_intercept', 'ANPP_intercept', 'overall_intercept'),
                   labels=c('Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall', 'Gamma', 'ANPP', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=3.5), linetype='dashed') +
  geom_vline(aes(xintercept=6.5), linetype='dashed') +
  coord_flip() +
  ggtitle('Richness Difference') +
  annotate('text', x=9.2, y=-0.8, label='(b)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 1600x1000




###by magnitude of resource manipulated---------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip, plot_mani)%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

rawTrt <- rawData%>%
  left_join(trtDetail)

#N addition (Figure 4)
Nmean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\n_mean_posteriors.csv', comment.char='#')
NmeanMean <- as.data.frame(colMeans(Nmean))%>%
  add_rownames('parameter')
names(NmeanMean)[names(NmeanMean) == 'colMeans(Nmean)'] <- 'mean'
NmeanSD <- as.data.frame(colSd(Nmean))%>%
  add_rownames('parameter')
names(NmeanSD)[names(NmeanSD) == 'colSd(Nmean)'] <- 'sd'
NmeanOverall <- NmeanMean%>%
  left_join(NmeanSD)

meanNPlotFinal <- ggplot(data=subset(rawTrt, n>0&plot_mani==1), aes(x=n, y=mean_change)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(name='Overall Community Difference') +
  stat_function(fun=function(x){(0.2003873 + 0.003866423*x)}, size=5) +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.7, label='(c)', size=12, hjust='left')

Nrichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\n_richness_posteriors.csv', comment.char='#')
NrichnessMean <- as.data.frame(colMeans(Nrichness))%>%
  add_rownames('parameter')
names(NrichnessMean)[names(NrichnessMean) == 'colMeans(Nrichness)'] <- 'mean'
NrichnessSD <- as.data.frame(colSd(Nrichness))%>%
  add_rownames('parameter')
names(NrichnessSD)[names(NrichnessSD) == 'colSd(Nrichness)'] <- 'sd'
NrichnessOverall <- NrichnessMean%>%
  left_join(NrichnessSD)

richnessNPlotFinal <- ggplot(data=subset(rawTrt, n>0&plot_mani==1), aes(x=n, y=S_PC, color=MAP)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(name='Richness Difference') +
  stat_function(fun=function(x){(0.04821532 - 0.0001258736*1000 - 0.01420553*x + 0.00001198415*1000*x)}, size=5, color='#4793CF') +
  stat_function(fun=function(x){(0.04821532 - 0.0001258736*600 - 0.01420553*x + 0.00001198415*600*x)}, size=5, color='#2D5E88') +
  stat_function(fun=function(x){(0.04821532 - 0.0001258736*200 - 0.01420553*x + 0.00001198415*200*x)}, size=5, color='#153049') +
  xlab('') +
  annotate('text', x=0.4, y=0.7, label='(a)', size=12, hjust='left') +
  theme(legend.position=c(0.02,0.02), legend.justification=c(0,0))


#H2O change (Figure S3)
H2Omean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\h20_mean_posteriors.csv', comment.char='#')
H2OmeanMean <- as.data.frame(colMeans(H2Omean))%>%
  add_rownames('parameter')
names(H2OmeanMean)[names(H2OmeanMean) == 'colMeans(H2Omean)'] <- 'mean'
H2OmeanSD <- as.data.frame(colSd(H2Omean))%>%
  add_rownames('parameter')
names(H2OmeanSD)[names(H2OmeanSD) == 'colSd(H2Omean)'] <- 'sd'
H2OmeanOverall <- H2OmeanMean%>%
  left_join(H2OmeanSD)

meanH2OPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=mean_change)) +
  geom_point(size=5) +
  scale_y_continuous(name='') +
  stat_function(fun=function(x){(0.1820251 + 0.0002544999*x)}, size=5) +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.65, label='(d)', size=12, hjust='left')

H2Orichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\h20_richness_posteriors.csv', comment.char='#')
H2OrichnessMean <- as.data.frame(colMeans(H2Orichness))%>%
  add_rownames('parameter')
names(H2OrichnessMean)[names(H2OrichnessMean) == 'colMeans(H2Orichness)'] <- 'mean'
H2OrichnessSD <- as.data.frame(colSd(H2Orichness))%>%
  add_rownames('parameter')
names(H2OrichnessSD)[names(H2OrichnessSD) == 'colSd(H2Orichness)'] <- 'sd'
H2OrichnessOverall <- H2OrichnessMean%>%
  left_join(H2OrichnessSD)

richnessH2OPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=S_PC)) +
  geom_point(size=5) +
  scale_y_continuous(name='') +
  xlab('') +
  annotate('text', x=-80, y=0.6, label='(b)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(richnessH2OPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanH2OPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600



