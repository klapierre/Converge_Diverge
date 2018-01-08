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

rawData <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\for resource type analysis_nov2017\\ForAnalysis_allAnalysis.csv')

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
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_3.csv', comment.char='#')
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
# write.csv(chainsCommunity2, 'bayesian_output_summary_expinteractions_12182017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteractions_12182017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,3470:5941]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 3470:5941])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,3470:5941]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 3470:5941])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_12182017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_12182017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=8, 7, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1721255)+(0.3228545), (intercept+linear*7+quadratic*7^2)*(0.2468055)+(-0.07086403)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1721255)+(0.3228545),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2468055)+(-0.07086403)))%>%
  mutate(color=ifelse(trt_type=='CO2', '#fabebe44', ifelse(trt_type=='N', '#3cb44b44', ifelse(trt_type=='P', '#aa6e2844', ifelse(trt_type=='drought', '#ffe11944', ifelse(trt_type=='irr', '#0082c844', ifelse(trt_type=='precip_vari', '#46f0f044', ifelse(trt_type=='burn', '#f5823144', ifelse(trt_type=='mow_clip', '#00808044', ifelse(trt_type=='herb_rem', '#f032e644', ifelse(trt_type=='temp', '#e6194b44', ifelse(trt_type=='plant_mani', '#911e6444', ifelse(trt_type=='R*burn', '#f5823144', ifelse(trt_type=='R*mow_clip', '#00808044', ifelse(trt_type=='R*herb_rem', '#f032e644', ifelse(trt_type=='R*temp', '#e6194b44', ifelse(trt_type=='R*plant_mani', '#911e6444', ifelse(trt_type=='R*R', '#80800044', ifelse(trt_type=='N*N', '#80000044', ifelse(trt_type=='all_resource', '#00008044', '#80808044'))))))))))))))))))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,',
                       '*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,'),
         curve5='), colour=grey)+',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_12182017.csv', row.names=F)

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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.69043731121695 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.487355353586635 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.32338948426585*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.817661631996 + 0.3566719400115*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(1.1191801886805 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.60694045551955 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.237183228925761*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.31608126157185*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.263094447150895*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7167174853536 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.57698644943935 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.97987667835 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.71646141761 + 0.16735275999705*x + -0.0094092267835572*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7349798372026 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.61851559716385 + 0*x + -0.022861438127795*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.23715365776265*x + -0.029169744475785*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.026024840399752*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.5176154270717 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.7776647094121 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.491072045873455 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.6174821159114 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.63570846909545 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.63872990959125 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.698807676315 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.30089578375 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.681113444447 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.42112245265 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.24538404405 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.2253270667 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.598838996975 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.1202431239 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.78350757719925 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.898537339735 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.09348761805 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.364868524717375*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.8543537807635 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.77791609556505 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.658822602224675 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.950044252795 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.875782409555 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-0.83801011934 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.10991954605 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.0168676507 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.25271127685 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.15264310185 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.109714347005 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.6983669429665 + 0.1974998293511*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.527388190279135 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.47437406376425 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.7685600348933 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.6642529722215 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.80750899658 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.10622133535 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.09046662375 + 0.15041580179235*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.851011971115 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.9150445331 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.766504146315 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.735468878153 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.5171589800247 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.7889815924606 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.840607535985 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.62470889431335 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0.024020210497*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.6830194833057 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.207803029565077*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.165785790419305*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.21071860560255*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.775613539262 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.71674478392505 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.8614745384 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.92419443326201 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.98941706035 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.0501798527 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.5967368564397 + 0.33567000335465*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.06185060505 + 0.25566868006513*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.66196697710635 + 0.3161992714019*x + -0.0184143013286425*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.7167408120394 + 0.2311407244806*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.23359893227735*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.53525462669445 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.48061364419634 + 0.3968239870127*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2768731852488*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.6048605411617 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.6422261666313 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.049381492275 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.10734473155 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.9333661486045 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.16030889755 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.0449288903 + 0.179426213174455*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.3705836852 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.0374454522 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.102278079988 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-0.880540244415 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-1.41579088195 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.5096354638015 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.814423212975 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.521109276687485 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7034625870535 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.4057240403 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.5770222841945 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.915455236635 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.749206837 + 0.229776653855*x + -0.012381044462105*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.00010388275 + 0.215083888496595*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.93924227115 + 0.22631096764015*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.94618765935 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.682209976728 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7338521210685 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.0516859207 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.54497728962595 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.5320322901295 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.6808061178615 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.7958799243615 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.4898075667104 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.711584493637 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.4917930003854 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.6446590948962 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.219661200709275*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.88622643497 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.5243387670915 + 0*x + -0.0112698077950582*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.87431307875 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.92355605332875 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.47153271735 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.039221194 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.60573919820215 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.88322881095 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.97231904355 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2393812790095*x + -0.0131359333014*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.6235595564845 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.8484690757515 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.84064324398 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.1319239574 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.910224661495 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.8061534457642 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.97151988045 + 0.3602064039235*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.9250994974 + 0.2473923312316*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.56738276454355 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7288472617913 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.7183500081972 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.56152921337215 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.74746558424625 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.3527849770835*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.990562208555 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.78563949991011 + 0.6903147863758*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.65365344942045 + 0*x + -0.046437623548*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.52950566044965 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2758689644043*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.66640993333655 + 0.25954644125476*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.7828037927042 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.9587822190913 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(1.2355192373 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.458026412662785*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.25796447413938*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.30894067070875*x + -0.019133741261905*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.72465967895355 + 0.30673971417*x + -0.01878378712957*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.91116377895 + 0.150681935276684*x + -0.010534619486693*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.304712253306595*x + -0.03791501560504*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.73826058964605 + 0.243708135856545*x + -0.03393144281192*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.3289679236531*x + -0.0390883411166115*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.3337610540528*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.66930414735715 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.7167838195755 + 0.2527114762505*x + -0.03352985405047*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.772250481102 + 0.235632576927235*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.4365376177 + 0*x + -0.0267835718576135*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.37768928955 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.4203504276 + 0.23837515579094*x + -0.0303564082255*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.654920840651 + 0.45333182108*x + -0.0295236286885914*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.4733828323428 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.68102404296605 + 0.26826984652902*x + -0.0336926793055*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.744235262491 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.757995249349 + 0.27167600849255*x + -0.0338848857426*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.5630632601351 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.8439293834665 + 0*x + -0.0315570570702*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.893254258225 + 0*x + -0.0251844945357677*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.0339021753481845*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.02646368515 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.21927366455 + 0.36320318885945*x + -0.0404023208682*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.2154616054 + 0.36781822208205*x + -0.03018086778779*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.0898894758 + 0.292609686298745*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.0641321468 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.6386264192755 + 0.266807291345*x + -0.01754289980136*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.737245221266 + 0.24700708597351*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.685152300385 + 0.2269437776914*x + -0.0128248475993295*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.614665102155982 + 0.19041409824649*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7991010125325 + 0.185103043155735*x + -0.0137741022138705*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.526673080765047 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.34348885951375*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.31869067613436*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.3356323239342*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.6697639490345 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.1769323259418*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2591970677275*x + -0.01196248842966*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.173988226476005*x + -0.00982432738109*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.25648659715*x + -0.01527237951*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.1895323234756*x + -0.0087690408182625*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.152147935903164*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.716148452074 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.5400329729868 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.10016523425 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.27567527631245*x + -0.02203368332045*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2086908346772*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2380006328869*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.93221768075 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0.0201364228890232*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0.0236181540429365*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.318647652995*x + -0.01553273230491*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.30685676837*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.69513964006696 + 0.381546693887*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.32587207425 + 0.3996713427*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.3614243472 + 0.444275689345*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.0990596336 + 0.32130497144*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.21019242305 + 0.37429493845*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.64460397632745 + 0.29524866359*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.51637103129005 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.5617800630963 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.80896929391 + 0.40051685875*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.17026212025 + 0.2123326365242*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.1477274107 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.97859253705 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.96527264164 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.7904328213695 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.892275026765 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.52375453435455 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.530711831233385 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.36987713575 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-1.07502185425 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.913489397215 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.63231512668725 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.36055837232*x + -0.026194858790555*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0.48101691025125*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-0.535078572941465 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.5195728602405 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.59457214514215 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-0.535884090366335 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.3695322817 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.39389052155 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(-1.48520524015 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.5387708546 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.47751306285 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.39383567115 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.5215790461035 + 0.151229805695*x + -0.00257351270136595*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.77637080825 + 0.211121679745*x + -0.0048093035902*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.88639545415 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.3436785522*x + -0.0082938813925*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.4165367187*x + -0.010034934042*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.6162424979283 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(1.1440016633 + 0.087238431329945*x + -0.0032332993285645*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.9599059536 + 0.102933831433*x + -0.00221609367488985*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0.07541683194112*x + -0.0037322329220365*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.81246598925 + 0.075435583544738*x + -0.00436186292571325*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.759435094355 + 0.230029017165*x + -0.0099035806545*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.784857823804 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.6254076601515 + 0.29955972335*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.267536370071112*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.6386416335*x + -0.0305260782015*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.94152752705 + 0.6021986922*x + -0.033455397257*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.7000624158858 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.4561198796096 + 0.37665940585*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0.49133040785*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.852938524945 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-1.15231360015 + 0.2937838359532*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0.3326399384625*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.444124129267 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.5815827455098 + 0.196299432913375*x + -0.0128256281031585*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2232576613375*x + -0.0159387374779275*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-0.5073815823362 + 0.202416491471565*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.1864482984953*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2746712250902*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.27301533798643*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-0.44623731861205 + 0.391250371175*x + -0.021008087372965*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.34569170134*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.20476090186308*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.221563636917205*x + -0.0069654698212515*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2231590494965*x + -0.00885700207324*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-0.494792465008615 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-0.88153179444775 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.118530746395 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.34345326406225*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.6520963868542 + 0.269620746128*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7569166499425 + 0.3955000513505*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.95152895186 + 0.298790664699*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0.235502724058202*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0.213324278925905*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0.311250035719*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2779134641105*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.53106350120295*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.53292923438016*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.0090962464485 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.75588285011625 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.537260957165095 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.7096931637458 + 0.51942773866837*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.4803580095765*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.628612064691 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,2), colour='grey')+
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.36334550 + 0.18923200*x + -0.00000000*x^2)*(0.1721255)+(0.3228545)}, size=5, xlim=c(0,7), colour='black')

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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.57625688691969 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.33942822818245*x + 0.032863792120425*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.4537920097225*x + 0.040735193075*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.13865660675241*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.14738085955005*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.7283473971925 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.2852149667119*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.7108334124193 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0.63920451024728 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.265311112014145*x + 0.024921058888925*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(1.09935510354 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.5391192921812 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.710188458237 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.76305150652445 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.85444182040835 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0.76070867262805 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0.60686992553505 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0.6256202598191 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0.9483809161315 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + -0.241216098315005*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.746374829685 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.7678721258584 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.926107752578 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.9558352019855 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.74784332350615 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.57810259550805 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.56359778978449 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.59648471437529 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.53233845563225 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.6379021215006 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.560976652132735 + -0.3229174188313*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.820899702 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.05727087374635 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.61020239495*x + 0.0462573622661376*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.9501890947975 + -0.477027663109*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.565392490215*x + 0.04518641353115*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.13414573764905 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + -0.457351282390625*x + 0.047259460447935*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.7316280220269 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0.70906296795035 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.9306652633307 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.55644098641925 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.18749626131687*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.981636141735 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.62583454943268 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.217283882190345*x + 0.015858578795875*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.8270878351735 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.6828063644753 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.53737194890094 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.53984486138955 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.16036083169 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(-1.009744983859 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.8821125977145 + 0.4726732428289*x + -0.0473246354625*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.0299682145046495*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.87221626985565 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.6783476199645 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.83649251511775 + 0*x + -0.023098288428885*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.16708177501263*x + 0.010806837892591*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.7606806867328 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.5918048425367 + 0.3207388350813*x + -0.037021397447405*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.706188021221325 + 0.290334417429026*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.67314653072992 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.563429330983625 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.724987243199225 + 0*x + 0.0535308418931*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(1.1850108471 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.90828387311325 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.7055504861203 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.591876041872405 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.9570846267125 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.57032153928744 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.3496769654416*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.358275058261075*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0.8433001444538*x + -0.24836936172*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.279576824439615*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.5315913841617 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.3473550947 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.21112745795 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.19772256195 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(1.2519671128 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.9686668372995 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.8086186644113 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0.54347044818215 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0.7610618812816 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.017371140499153*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.0156797939760404*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.645791218235625 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.0293888546349585*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.530863950249275 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.2287781201542*x + 0.015551438206*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.231155838593555*x + 0.0135046325849*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.841851994079 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.738141645348725 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.3079267425748*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.2962490954865*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.240225951815272*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.5915439428743 + 0.31893580110975*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0.320281806681154*x + -0.036660408025*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0.51486427175465*x + -0.05449685022215*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0.5628897051263*x + -0.06008008002915*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(1.14171476075 + 0.52986093248335*x + -0.06799207075*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.8813409673082 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0.61886367895*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(1.2438132054 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.994387101882 + 0*x + 0.6012216161*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.829288758995 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.5966191021085 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.5774563714655 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.6396906394955 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.817865165614 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.600493805602695 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.5863379205353 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0.825979513489 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0.663311387794 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(1.59661965125 + 0.8149584306498*x + -0.080607672481045*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0.8213640474935*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0.5079846106642 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.075048315466802*x + 0.0027480275506404*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.30158857175 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-1.80797373715 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.989658617725 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.57602782264545 + -0.32103210475*x + 0.008366788679*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.25441766875*x + 0.00626624766885*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(-0.90634998928 + -0.141825403862*x + 0.0043491682147917*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,5), colour='grey')+
  stat_function(fun=function(x){(-1.1225760543 + -0.12573646370712*x + 0.0033267611004675*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,3), colour='grey')+
  stat_function(fun=function(x){(0 + -0.177653482091*x + 0.00480463625777065*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.8217690681*x + 0.06433949108*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(-0.7069759015485 + -0.8239339929*x + 0.06381660758*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.33054722569603*x + 0.034709732901273*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.686280686*x + 0.056841275925*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + -0.8163433263*x + 0.05475926224*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,7), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + -0.308513175538*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(-0.71874584951185 + -0.458190074375*x + 0.030304547424365*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + -0.2976977759841*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + -0.3254276296466*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + -0.3161664283767*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + -0.3261898121265*x + 0.025299830414905*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0.0113142221717199*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + -0.24419974050475*x + 0.01248881014917*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + -0.18802896651275*x + 0.0127287916949605*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + -0.4207471531*x + 0.02648606954*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0.73851035600265 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0.609921640813945 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0.81167897294*x + -0.075600373333575*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,4), colour='grey')+
  stat_function(fun=function(x){(1.36949989255 + 0.7270234762*x + -0.10386264202*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(1.13771274205 + 0.28810688626635*x + -0.04421706428345*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.8862194924995 + 0.6438528377*x + -0.095700136335*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + -0.036932028451102*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,6), colour='grey')+
  stat_function(fun=function(x){(0 + -0.27827809467936*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.573551096368165 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0.9797461640375 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0.58575624165175*x + -0.119312336085*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,2), colour='grey')+
  #overall line
  stat_function(fun=function(x){(0.28671150 + 0*x + 0*x^2)*0.2468055 + -0.07086403}, size=5, xlim=c(0,7), colour='black')

# print(richnessPlot) #export at 1200x1000

#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,1)))
print(richnessPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(meanPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#export at 1200 x 2400



# ##summary stats from bayesian output --------------------------------------------------------
# # gather summary stats needed and relabel them
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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_12182017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_12182017.csv')

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
  ggtitle('Community\nDifference') +
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
  ggtitle('Richness\nDifference') +
  annotate('text', x=9.2, y=-0.8, label='(b)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 1600x1000