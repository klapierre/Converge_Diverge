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
setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

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

#laptop
rawData <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10 yr\\ForAnalysis_allAnalysis10yr.csv')

test <- rawData%>%
  select(site_code, project_name, community_type, treatment, experiment_length)%>%
  unique()

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
# chains1 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10 yr\\simple_mod_effs_coding_interactions_10yr_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10 yr\\simple_mod_effs_coding_interactions_10yr_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10 yr\\simple_mod_effs_coding_interactions_10yr_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\10 yr\\simple_mod_effs_coding_interactions_10yr_3.csv', comment.char='#')
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
# write.csv(chainsCommunity2, 'bayesian_output_summary_expinteractions_10yr_01122018.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteractions_10yr_01122018.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,3326:5683]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 3326:5683])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,3326:5683]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 3326:5683])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_10yr_01122018.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_10yr_01122018.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=10, 9, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1510305)+(0.2985893), (intercept+linear*7+quadratic*7^2)*(0.233634)+(-0.04097597)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1510305)+(0.2985893),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.233634)+(-0.04097597)))%>%
  mutate(color=ifelse(anpp<281&variable=='mean', '#1104DC44', ifelse(anpp<561&anpp>280&variable=='mean', '#4403AE55', ifelse(anpp<841&anpp>560&variable=='mean', '#77038166', ifelse(anpp>840&variable=='mean', '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,',
                       '*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,'),
         curve5='), colour=',
         curve6=') +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_10yr_01172018.csv', row.names=F)

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
  scale_x_continuous(limits=c(0,9), breaks=seq(0,9,1), labels=seq(1,10,1)) +
  ylim(-10,10) +
xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  annotate('text', x=0, y=1, label='(b)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.60533764232735 + 0.39440937129465*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.827115739713 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0.73326946949225 + 0.31851865188571*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(1.093665306894 + 0.452045076072685*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(1.295619442275 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0.58849243317897 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0.7895547922595 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.35861277913223*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.4995759667905*x + -0.0434440472436645*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.49878201240375*x + -0.047110248281135*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7371223907335 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6436487017535 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.980478875465 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7442373440372 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.738351874595 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.527852323948075 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.50787768379212 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.61545214376 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.04606731329546*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.353627637020995*x + -0.048430833981572*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.498539449056743 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6463051213236 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.364491074983435*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.812055616583 + 0.3954719775*x + -0.044267645006*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.29393120515 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.783890403796935 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.3694631383 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.24693025555 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.25294252975 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4256746923968*x + -0.0476555635484*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.50716374312345 + 0.34108150343915*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.874036892245 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.033989360928 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.7753349153625 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.8657155310475 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.0536418464245 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.797592130931 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.7122770264479 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.643405513731 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.84618349803 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.913501673246 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(-0.9389835501575 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(-1.18847137865 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(-1.08586393225 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.2917284031 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.14361531235 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.10844279607 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.758934161204 + 0.32486217252785*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7232079756649 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.896701391265 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.878705020395 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.00719110247 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8830252773515 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.930353471835 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.925160256215 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.867472160965 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.65005645279785 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7095850210316 + 0.405413870568*x + -0.032603491303008*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9054706816 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.73459101617945 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7137522118872 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.03280267904982*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.801595288924 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.753189046355 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.85831404434 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.882234641491 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.93134782278 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0611244674 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.55929636895525 + 0.394431349454045*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.19615933285 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.866336064795 + 0.4948111734965*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.960399870455 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0.820375057433 + 0.4395256126858*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.649961792085845 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(-0.80629969696255 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(-1.23985259255 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.015325927875 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.274980725 + 0.33959900752465*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0595984952 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.18459243525 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.08341324605 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.9607163888 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.95141445137 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.49417604467611 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.3305887605 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.46619625079959 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6993683377025 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.69947729389085 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.081617210926715*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.5947236649459 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.2874414802 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.788955191045 + 0.372620986475*x + -0.0501593744453*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.9824996878 + 0.28564624365952*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.0099569314 + 0.297105058547297*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.9259050893 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7951118705815 + 0.385848178871145*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.74842483044995 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.5411226465305*x + -0.049340780791145*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.01441457885 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6967196732455 + 0.2729843707955*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.567471360086 + 0.35248551143712*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.426787280084135 + 0.26995470709448*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.35275928465854*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.853684578915 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.90950886232 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6503428621265 + 0.3087188122313*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.53107270587485 + 0.26474103239861*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.63108990339075 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.561747633946045 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.784478504077 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.42105864418255*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.0236695622985 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.3748650992054*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.66166427256485 + 0.34847774026515*x + -0.03458466858257*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8699644485145 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7533338426805 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.3608359276 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.98784210272 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.4799135430978 + 0.23421535261179*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.320140507516515*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.45609824775445*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.87832033909 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.81959873653 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.58273332781425 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.59878825692545 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.871853472775 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.1137465009 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.850467698422 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6969027094793 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.863879914384 + 0.29459968450132*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.8643891165385 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.370322655030955*x + -0.05571828663166*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.751077114064 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.80674092167 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.4785279413953 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6602808250875 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4423324836225*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7343330371368 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.283847642013402*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0955765909391*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.005386461535 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.59570320657135 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.53316175295165 + 0.436127280891745*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.8877428438895 + 0.4539659732905*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(1.6332761351 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(1.1877617791795 + 0.83083224684654*x + -0.18470456139*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(1.63360161985 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0.64761194005615 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0.5495358029754 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0.633836217946345 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.5085468745486*x + -0.044659643566975*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0.57001517707065 + 0.4195857088859*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.496450173332*x + -0.048861931000195*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6270941606125 + 0.3276318506972*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.82848832877 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6049557479652 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.32545337691437*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.87760331323 + 0.413785408078*x + -0.05291562290835*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7835218071225 + 0.29730342892731*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.892674434535 + 0.480553662715*x + -0.051702355710315*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.41752066955 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.46818609875 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.3483751668 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.5474717288152 + 0.357284219571515*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.561938910332785 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.61684896354465 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.65599918561013 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.5166097855992 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8240213119201 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.803302565763 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.02382982705 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.0793548651 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.22444495161 + 0.51188201406485*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.79126400831 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.876647832712 + 0.514799564895*x + -0.0434495686702*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.51483797496935 + 0.30879261166655*x + -0.0261078910231905*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.51200907769305 + 0.221954338670145*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.53384414690585 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.560810144706*x + -0.04432814562089*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.518206713656*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.48672263993535*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.5262251830935*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6836648678335 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.76483860959 + 0.360956198735703*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6769398391651 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,8), colour='#DD032688') +
  stat_function(fun=function(x){(-1.03140347685 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.434851160194*x + -0.03721002615624*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4401203884365*x + -0.046944190718265*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.566331452079 + 0.58520259085*x + -0.0509792509775*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.954498480635 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + -0.28833447473531*x + 0.041261799610195*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.46189819248*x + -0.027143990757249*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.61115384711675 + 0.3875674537213*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.16319620705 + 0.268912093708082*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.18225224665 + 0.26342584588163*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.23863672115 + 0.455677677725*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.2099008303 + 0.2987868458638*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0.37912194922112*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + -0.97327076766555*x + 0.39789193370885*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0.48896937409765 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.842262194535 + 0.29915527595295*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.94117088955 + 0.52220151453*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.99047917365 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.09053810755 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.99124466685 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.870572594386 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(-0.7184345158215 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(-0.51439439012415 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(-0.45340181964875 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#77038166') +
  stat_function(fun=function(x){(-0.53064447546495 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.4440690316 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.95523686948 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.905147227765 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.5756412898083 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.3587282396582*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0.715945223336 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.3569759481375*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.70079405528965 + 0.32988025442655*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.5285901896845 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.45277089225 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.7709888889 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.52041366665 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.7464916166 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.61206088375 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.6376554332 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.675189376634395 + 0.313464895295275*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.521726575835*x + -0.030451523014467*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.61142007564*x + -0.04387684138841*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.79570815663985 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7670407229435 + 0.49090488116*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7727884635575 + 0.6408671627*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7760677642625 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.87104111655*x + -0.075040952425*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.661464401*x + -0.04453953857294*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.30076656352006*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.647921437515*x + -0.04253145363555*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.7853551365*x + -0.0572563277495*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(-0.56791334667335 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(-0.72048251865295 + 0.461774606517*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.27395770390915*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.84182326785*x + -0.0481481543868*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0.7334093495275 + 1.00327805615*x + -0.073716899665*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.5457170745952 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.672494156695 + 0.5783379434*x + -0.0287873323513935*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.75201486715*x + -0.0402708129537*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.1665069854 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.042886805375 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.816916053157 + 0.4989634432244*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + -0.4277223301138*x + 0.0786136802122*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.45084628401645 + 0.3502344362943*x + -0.033704443698755*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.47493662568561 + 0.237966274034347*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.383441352671*x + -0.03856959054325*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.272597932271005*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4222569985185*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4577884682005*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.595700674835*x + -0.04290060670995*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.280869443732082*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.68109541433815*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.37647118453185*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.5588877663009 + 0.319600524082243*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.870611217353 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.08159764765 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.61593614366*x + -0.037432467683106*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6828259932925 + 0.355694111841*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8698818388005 + 0.5070914396025*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.9575351702 + 0.351372801502*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.392371949037*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.44665165588825*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.62956781598925*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.61953127698415*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.8747907719775 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.6703064485844 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0.630972620727135 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.59679211782385 + 0.5501949698293*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.53776994910535*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1510305)+(0.2985893)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.5017765 + 0.2313645*x + -0.00000000*x^2)*(0.1510305)+(0.2985893)}, size=5, xlim=c(0,9), colour='black')

# print(meanPlot) #export at 1200x1000



#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-1.0,1.0))  +
  scale_x_continuous(limits=c(0,9), breaks=seq(0,9,1), labels=seq(1,10,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('') +
  ylab('Richness Difference') +
  annotate('text', x=0, y=1.0, label='(a)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.50702924653215*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.685714530765*x + 0.05598079926714*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.45665996068675*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.62200443069425 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.44004572735759*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.7280192103535 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.965270853071 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.76308097634765 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.884885056219 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.8338865465215 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.36143144625616*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.7923551686825 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.9202082549045 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.9279665199715 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.598862880876555 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.64193439492991 + -0.502917322593286*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.705859918834635 + -0.4959612267041*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.95211110795 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.82855407345955 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6387141666309*x + 0.057422526853265*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.73679520704*x + 0.06484240429463*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.8626256691345 + -0.4963906602433*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5769799493975*x + 0.06751507120766*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.050978310271 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.67946873643185 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.63157004520635 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.9558830500555 + -0.57303197023269*x + 0.086421578301745*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6000267432378 + 0.4515278974525*x + -0.053205777349*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6842737812552 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5961449621039 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0331453135722806*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4475312488475*x + -0.05122594748755*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.5556689299 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.804071959463 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.07663910664485 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.34680827765 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5247214118171*x + -0.08214444034248*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.8760818238935 + 0.5788186976175*x + -0.063378305280335*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.89022224378263*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.328933218031825*x + -0.05186442548215*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5697673148955*x + -0.0732967786302*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5161177926656*x + -0.066567996365625*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04453110391053*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0376139576724825*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.547182978142*x + -0.0636626106545*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.72348223365*x + -0.085933461415*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05031301566939*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.71642787846825 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.58865257536095 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.411768163618935*x + 0.06768316963705*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.906821140277 + 0.411955212571315*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.8199638112108 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.19924180686 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6851114968515 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3259052957444*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.471849129183635*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.60699537495775 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.2757896201 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.093448424145 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.20651679581 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.08946628772365 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.9282446975825 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.9617337571185 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.73192418007745*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36084469860002*x + -0.038165229792905*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0365966496039585*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5149447595524*x + -0.055908718627735*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5449377555025*x + -0.0695494902373858*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.572352364314535*x + -0.0698942595368*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.58056728110275*x + -0.06546432102285*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.61870133064005 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7905559980125 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.582850352403*x + 0.043128050729565*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.41521357379545*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0464355021333*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0567563002286865*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.45237269271425*x + -0.050288508107626*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.58618197529*x + -0.0661768475575*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6018627772765*x + -0.061720714046*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.65457369623*x + -0.0698009911365*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.97424030133 + 0.7436089928715*x + -0.1013650628265*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9217337110657 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.85844384865 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.24948061016 + 0.7011075609296*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.57127664253265*x + 0.0486295809016755*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.51973435101255*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.469546046273555*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.863449679304 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.5659698794845 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.74640511048585 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.61887880065985 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.51772162498785*x + -0.0711981963287155*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.43391722159 + 1.3097076532*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.998719067013*x + -0.125279062264*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04618019729299*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.98001438745*x + 0.071035575431*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -1.1743341724*x + 0.093240227175*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7937550585405 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.98802860424955 + -0.64132765813465*x + 0.0341150969835845*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.81578444275*x + 0.05708260977603*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.7943060466*x + 0.06636128239035*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.66621769425*x + 0.0504673974166*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.63100007922415 + -0.562848487455*x + 0.0526273160337*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.631657871925*x + 0.061000049101*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -1.01403498515*x + 0.0754607486215*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7894628432734 + -0.9682647018*x + 0.080097321845*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.89846981005*x + 0.07541383032*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -1.07400916515*x + 0.087761913065*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5203461522745*x + 0.03895368668166*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.64002033168224 + -0.640214832565*x + 0.04848046208715*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.425424565323465*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05179857766092*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04352495492174*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05754585593036*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.082158526321713*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03574219425647*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.290329157268525*x + -0.0341942094761765*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.28399903499745*x + -0.039769767700505*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.82790654658895 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.878982389495*x + -0.08321555626888*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(1.113766647035 + 0.99408741675*x + -0.13743172293*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9641886199635 + 0.43233863650615*x + -0.061594321813*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.65150329637813 + 0.8557974765*x + -0.122613081605*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32666550960275*x + -0.048017555065645*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.377511084849345*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4607927154788*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.627056782238625 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.038735140334 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.233634)+(-0.04097597)}, size=2, xlim=c(0,2), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.2279295 + 0*x + 0*x^2)*0.233634 + -0.04097597}, size=5, xlim=c(0,9), colour='black')

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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_10yr_01122018.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_10yr_01122018.csv')

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
  annotate('text', x=9.2, y=-0.8, label='(b)', size=10, hjust='left')

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
  annotate('text', x=9.2, y=-0.8, label='(a)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
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



