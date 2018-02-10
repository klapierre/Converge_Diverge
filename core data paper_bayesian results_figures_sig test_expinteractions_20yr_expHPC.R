library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)
library(nlme)

###note, for this analysis S_PC is really expH_PC

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

rawData <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\20yr_expH_PC\\ForAnalysis_allAnalysis20yr_expH_PC.csv')

test <- rawData%>%
  select(site_code, project_name, community_type, treatment, experiment_length)%>%
  unique()

rawData2<- rawData%>%
  left_join(trtInfo)%>%
  filter(plot_mani<6, anpp!='NA', treatment_year!=0)%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_rich=mean(S_PC), std_rich=sd(S_PC)) #to backtransform

###prelim figs - test
ggplot(rawData, aes(x=treatment_year, y=S_PC)) +
  geom_point() +
  geom_smooth(method='lm', se=F)

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
# chains1 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\20yr_expH_PC\\simple_mod_effs_coding_interactions_20yr_expH_PC_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\20yr_expH_PC\\simple_mod_effs_coding_interactions_20yr_expH_PC_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\20yr_expH_PC\\simple_mod_effs_coding_interactions_20yr_expH_PC_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\La Pierre_comm difference_final model results_01122018\\20yr_expH_PC\\simple_mod_effs_coding_interactions_20yr_expH_PC_3.csv', comment.char='#')
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
# write.csv(chainsCommunity2, 'bayesian_output_summary_expinteractions_20yr_expH_PC_01122018.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_expinteractions_20yr_expH_PC_01122018.csv')

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
# write.csv(chainsFinal, 'bayesian_output_mean sd_expinteractions_20yr_expH_PC_01122018.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_expinteractions_20yr_expH_PC_01122018.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=20, 19, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.165778)+(0.3160571), (intercept+linear*7+quadratic*7^2)*(0.4724521)+(-0.1514531)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.165778)+(0.3160571),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.4724521)+(-0.1514531)))%>%
  mutate(color=ifelse(variable=='mean'&rrich<31, '#1104DC44', ifelse(variable=='mean'&rrich<51&rrich>30, '#4403AE55', ifelse(variable=='mean'&rrich<71&rrich>50, '#77038166', ifelse(variable=='mean'&rrich>70, '#DD032688', 'grey')))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,',
                       '*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,'),
         curve5='), colour=',
         curve6=') +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_expinteractions_20yr_expH_PC_01172018.csv', row.names=F)

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
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  ylim(-10,10) +
xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  annotate('text', x=0, y=1, label='(b)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.489695734110546 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7135350845633 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.537731089069021 + 0.3063560521093*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.927572202845 + 0.3826344774115*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(1.176757560965 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.705130432811 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.27882671238948*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.37979809896*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.3483234031193*x + -0.020701138999906*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7561126308235 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.632400734517 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0260302004 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7369947346015 + 0.18778931241365*x + -0.0106680498661215*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.76958257341 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6250083717626 + 0*x + -0.029831626518589*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.3393865725661*x + -0.0418526151282*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.250686278577985*x + -0.0362897603014*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.55936711654154 + 0.18487356373109*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7186679621345 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(-0.5930801359121 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(-0.622223347686875 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.5867932385702 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.776148906055 + 0.178749700552*x + -0.0133830330226895*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-1.36715363705 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6979039100072 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-1.4262431278 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.25727342555 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.2004274142 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.2410024463816*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.446929577406758 + 0.215259172158896*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.583747861300475 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.085437442045 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.76762002920527 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.90350470011 + 0.3985922091496*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.100358883695 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.48126593612355*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.8365188623786 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.75508596480625 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.68119026690225 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.900327580222 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.837681349005 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(-0.92402136037 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(-1.2345849991 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#77038166') +
  stat_function(fun=function(x){(-1.12388525063 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,15), colour='#77038166') +
  stat_function(fun=function(x){(-1.254123825 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,15), colour='#77038166') +
  stat_function(fun=function(x){(-1.1250453875 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,15), colour='#77038166') +
  stat_function(fun=function(x){(-1.0963206319 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,15), colour='#77038166') +
  stat_function(fun=function(x){(-0.70735547496275 + 0.216425255828175*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,15), colour='#77038166') +
  stat_function(fun=function(x){(-0.547500951095854 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7288143590075 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6223095820013 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7310903329932 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0602405692 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.10581299125 + 0.1565944455825*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.88184968325 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.92684913422 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7258243864225 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.837603330105 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.55025002758829 + 0.203478794634*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.856224193865 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.76370533611 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.50211095354285 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0277950964615*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.57405444712391 + 0*x + 0.01801856364881*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,13), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.21806201517975*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,13), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,13), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.20519231171165*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.210919201682175*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(-0.7059520546315 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(-0.735688483419 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(-0.8500113178325 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(-0.98972059054 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(-1.0698677147 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(-1.11072131083 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.61600994501902 + 0.37104327645213*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-1.04460995365 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6141297267762 + 0.32748027634866*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(-0.7137692938196 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.2454592511308*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#77038166') +
  stat_function(fun=function(x){(0.53295891129445 + 0.45524679138*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.24662722504945*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.608631177540735 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.61572596073645 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.129311223675 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.21433372805 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.9237414819 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.17343073115 + 0.19767282652168*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.1079812969 + 0.19641294854681*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.36855827065 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0776867076 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.12735837035 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.925099450695 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.29952242985 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.4850127798625 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.786890197468 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.60562392688269 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(-0.571731492717 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.76604514886495 + 0.22362611659297*x + -0.020591732464515*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.34487883175 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.5921783346571 + 0*x + -0.0229596369646*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.93732048785 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.18858630735789*x + -0.024324796752335*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.7816998377 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.49806471325875 + 0*x + -0.0186103670874875*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.05310003665 + 0.23110323995034*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9963334425 + 0.304301563995*x + -0.020873998494106*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.98153438365 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7665572895395 + 0.28068094996411*x + -0.02289149674226*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7905142754675 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.33849552520125*x + -0.03189307392545*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.05012624935 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6236366629925 + 0.1998080815323*x + -0.021048093923142*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.58895442597475 + 0.2248383667062*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6585685717743 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.754202131703 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.704352185877 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.50475669458962 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.53580897241583 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.234715972795005*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.8423863233925 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.5530919394078 + 0.1669765411181*x + -0.013390323879418*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9228508947 + 0.21427573349922*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9664050275 + 0.24806946477255*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.4459848946 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.027135867205 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6444077122105 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.30381945207925*x + -0.023195072161665*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.8480671796 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.01222181265 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.399143685145*x + -0.025000161054*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.19607472045014*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.625065795726 + 0*x + -0.0216704964501526*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,12), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.781577471866 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.8100240088245 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.984888687 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.9047346446915 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.7944427024236 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-1.02075196735 + 0.38868491564*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.86691092521605 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.27848585641857*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.5725173207375 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.72972707437335 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7235808411415 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.60955281678127 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.654593300048165 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.436328172236955*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.97407325012 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.729534323947 + 0.712252873389*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.577602303146585 + 0*x + -0.049301282796853*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,19), colour='#DD032688') +
  stat_function(fun=function(x){(-0.55289931296588 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.2879204286353*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(-0.62031960070015 + 0.24931170109354*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(0.929853440559 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#77038166') +
  stat_function(fun=function(x){(1.1172626581 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(1.34064223055 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.65074473886989 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,18), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,18), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.325949317096055*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.5945585621607*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.3300945231541*x + -0.025850223211265*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.3916558683733*x + -0.0283354698796*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6811752640935 + 0.316229375905*x + -0.019515509496*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,11), colour='#DD032688') +
  stat_function(fun=function(x){(-0.8259647763 + 0.139601492040485*x + -0.010233591029242*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.3591120775845*x + -0.04812873445535*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6563563507617 + 0.26643532620055*x + -0.040849824619415*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.4062370332455*x + -0.05020417422*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.390909193747*x + -0.029448483514305*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7008779300365 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.76503124567 + 0.2682554526998*x + -0.0353296173327*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7694801664955 + 0.258211583341*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-1.42817471485 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-1.392113254 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-1.3760031196 + 0*x + -0.0273214609440225*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6772027283826 + 0.461770002872*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.52996874771 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.6941177133605 + 0.3003533400129*x + -0.030556787701415*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.701745773175 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.7608664849335 + 0.294597491158185*x + -0.02949746577451*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.57231720031475 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.82516608362 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#DD032688') +
  stat_function(fun=function(x){(-0.8151120457645 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.92937945996 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.1988732333 + 0.33546155894095*x + -0.037070256673195*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-1.22389653865 + 0.363667969251965*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-1.1574680575 + 0.321933662066872*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(-1.029063478243 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6259043713049 + 0.287832841985*x + -0.019406929742111*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.76537184142185 + 0.2692042947284*x + -0.0162690399265703*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7189388522285 + 0.264564832675*x + -0.0159242449088835*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.61851579245515 + 0.1944027022942*x + -0.012623493073497*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#DD032688') +
  stat_function(fun=function(x){(-0.719240663747 + 0.15875991494976*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.315225147947475*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.313975793840365*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.3228221500309*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.291424402888615*x + -0.016864190878714*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.5543462976782 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.36621765863*x + -0.0206622924653*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0.83965895824 + 0.248538630667065*x + -0.01843047641995*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.47559754645*x + -0.032441253295*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.28185013215745*x + -0.016853038569765*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.46794887165*x + -0.0262134700875*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.2480842822923*x + -0.0104570654765168*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.650705007934 + 0.2573663261493*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0274196396 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.31369210012375*x + -0.025009517711106*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.202851146179491*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.25336145347073*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.94766616585 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0.02765760295098*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.35886078951*x + -0.0183113021322935*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.33091814612895*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6695831259247 + 0.357673800403*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(-1.28048523525 + 0.358774000892*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(-1.3067274846 + 0.3852655075245*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(-1.0716637597 + 0.31803217760905*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.25011009935 + 0.38638582362*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.668341666617 + 0.320785469978*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#DD032688') +
  stat_function(fun=function(x){(-0.4660092264228 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#DD032688') +
  stat_function(fun=function(x){(-0.905947162385 + 0.4638133093*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.2291684955 + 0.23042483567352*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.1692920683 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.00161796458 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.94473159055 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.834498491675 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(-0.8382465261305 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(-0.45082556662865 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#77038166') +
  stat_function(fun=function(x){(-1.47716058295 + 0*x + -0.0288523503654215*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#77038166') +
  stat_function(fun=function(x){(-1.14380918591 + 0.29690763076665*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#77038166') +
  stat_function(fun=function(x){(-0.932421497 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#77038166') +
  stat_function(fun=function(x){(-0.6654904575711 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.406097090515*x + -0.028060662095535*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,10), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.52935633941975*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,7), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.6745730991412 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.62229716630305 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.3890033659 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.49001964155 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.59905914205 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.6404775291 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.5818325122 + 0.2341049491105*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.33843957865 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.28073817191*x + -0.00946704268045*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.689260656755 + 0.281265965445*x + -0.007743980621*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.018731149635 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.0087434533 + 0.6662818383*x + -0.025709509655*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.7352105405515 + 0.6126127197*x + -0.0198984732155*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.83769935452 + 0.2439706982045*x + -0.0119386533783*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,5), colour='#DD032688') +
  stat_function(fun=function(x){(0.671712634502 + 0.272959673565*x + -0.012283237455*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.12681468948906*x + -0.0064753670593666*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.32352420295*x + -0.0188110190485*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,3), colour='#DD032688') +
  stat_function(fun=function(x){(0.5543940963166 + 0.334991715*x + -0.0153933618415*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.65375899549675 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.517162678628063 + 0.264637268722*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,9), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.2849913406755*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,16), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,16), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.75061309745*x + -0.0409273272715*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,16), colour='#1104DC44') +
  stat_function(fun=function(x){(0.862955234225 + 0.73579091635*x + -0.046616953595*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,8), colour='#77038166') +
  stat_function(fun=function(x){(-0.76940437559 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(-0.5557417657511 + 0.446288852*x + -0.018397191563151*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0.63536818175*x + -0.03150208957225*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(-1.05134845525 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(-1.218076719655 + 0.38018295608*x + -0.0202114181190663*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(-0.6344515763164 + 0.422535584836*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#77038166') +
  stat_function(fun=function(x){(0 + 0*x + 0.026461260197955*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0.028325574759675*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.491653545104585 + 0.182833170836395*x + -0.011967109653938*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.22807979081045*x + -0.016790443845655*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.4431344558207 + 0.198595971200645*x + -0.0120810137558614*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0.18821240085915*x + -0.013239905389675*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.3753199330858*x + -0.02662840281725*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.29764605816897*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.4212366981922 + 0.43658141765*x + -0.025505518627256*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.34821776976*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.27609427237395*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.234390504620499*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.360233149874*x + -0.017891577227828*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4266725948135*x + -0.024526225139145*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.2231430139854*x + -0.015525096301405*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.4139969455713*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.27901789687608*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.59273305501505 + 0.268664555798795*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.835203772811 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(-1.10311606825 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,4), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.395593436922*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.74278693612 + 0.3383657651281*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-0.7746473762495 + 0.40289085845645*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(-1.0358112244 + 0.349161096735*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#1104DC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.249198630261645*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,6), colour='#DD032688') +
  stat_function(fun=function(x){(0 + 0.346975642562*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.31076111808225*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.5502915571089*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.55634865090005*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.9439066029455 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.783145742032 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0.5362758322762 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(-0.6709054209667 + 0.52021871032197*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0.50472938718183*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.165778)+(0.3160571)}, size=2, xlim=c(0,2), colour='#4403AE55') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.38111500 + 0.20332350*x + -0.00000000*x^2)*(0.165778)+(0.3160571)}, size=5, xlim=c(0,19), colour='black')

# print(meanPlot) #export at 1200x1000



#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-1.0,2.0))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(4,19,5), labels=seq(5,20,5)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,1)) +
  xlab('') +
  ylab('Richness Difference') +
  annotate('text', x=0, y=2.0, label='(a)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.4506865112672 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4574040747555*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.7315970412*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.2439152402025*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6462171638149 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.62667267634005 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.33220281518715*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.24482685329805*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5734665594075 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.587888475059375 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1709706086951*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.44180263944185 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.44810944831665 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.4744254518973 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.42414131592706 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.46136789301782 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.26529134281195*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.57042074662805 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.46431108007575 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.09558578046748*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.089933233470435*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.091332940078482*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.08929058926599*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.617642025218645 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.45066017794565 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.49331069095008 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.58934896562545 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.021838618126886*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.55132295317125 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.55094023585055 + -0.235104683807905*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.423555082130876 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.482513068065604 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.5775758538387 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.53506303944087 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.44265404709895 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.46362113133191 + 0.20123943398691*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.4255508822085 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.45830556390164 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.49243558016325 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.52463614437815 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.08194587205 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.76147804180505 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.4423884571843 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.475972758187885 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.4581011225681 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17143485486494*x + 0.0120645400820165*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.4387726894534 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.425372234364635 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.55439314255 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.4801184562838 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.5745635842444 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.507295395962575 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.44780488272385 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.429792173264365 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.458376719558025 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.8622842599773 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.528158013200595 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.66672668900685 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.5630710816992 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.6089277261479 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.505791871967185 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.367794149675*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.503636586338385 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4755149116025*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.66710371475*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.89441334855*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.5825673516552 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + -0.47984349745*x + 0.031576270846*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.36980566090445*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.5190154761308 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.523220938010825 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4914371167752*x + 0.058321783126595*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.474174989104045 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.467193145345595 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.45416988335565 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.24804615663195*x + 0.016494763567705*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.22203089656682*x + 0.0146200284992705*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2578834176468*x + 0.0153703889068325*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.25299977311515*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.92591181975*x + 0.059013563885*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.43925812371325 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.644526148201 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.5305334832978 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.49734530788423 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6002127769405 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.54258200032885 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.489502564054185 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.60557132821175 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.58097305405375 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.589034103264 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.11787733780327*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.45032692333945 + 0.28914826234795*x + -0.02543705530767*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.533170095806205 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.4752315023088 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.4595902462914 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.46797255140585 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.65110324196095 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.751331915924835 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.54906531095335 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.940636340263 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.159188485728795*x + 0.0097910437710104*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.58406022483505 + -0.6209229082*x + 0.037540535775*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.04695077835 + -0.7307362013*x + 0.04301121575*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6352938313885 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.611295395366593 + -0.26191054486*x + 0.0098755282434*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.30588730107*x + 0.007099213722255*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6427519767988 + -0.4228483241*x + 0.022849717336*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.5780337621089 + -0.5580149172*x + 0.029339218645*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0085196929404415*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.00745785944138514*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.58507929305*x + 0.042720936929*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.81180855825*x + 0.058004214989*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.34790315013085*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5917855175*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.44759072865*x + 0.0389685716071*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.8880015683*x + 0.063440254155*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.297186502081695*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.46956196427285 + 0.35813302452305*x + -0.077288461675*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.59571001933521 + 0*x + -0.0608449478934*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.347656646425225*x + -0.07621645193*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.68328507046445 + 0*x + -0.075325805535*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.67853524879 + 0*x + -0.06008291356119*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.5130748457497 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.290777220916615*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.298891137447595*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.568685552515*x + 0.03451312526*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2766802705763*x + 0.0149974874962105*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.71291717435*x + 0.04416852316*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.730574320790055 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.43343193592727 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.72310407873365 + 0.35035220762233*x + -0.056115898495045*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5959688296335 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5670035613086 + 0*x + -0.0484694296717145*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.48140470949385 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5618067974985 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.528237341774725 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.659731486135 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.58430620037525 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.71414549097992 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.4724521)+(-0.1514531)}, size=2, xlim=c(0,2), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.29902950 + -0.12054050*x + 0*x^2)*0.4724521 + -0.1514531}, size=5, xlim=c(0,19), colour='black')

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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_expinteraction_20yr_expH_PC_01122018.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_expinteraction_20yr_expH_PC_01122018.csv')

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



