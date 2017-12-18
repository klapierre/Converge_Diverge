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
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich\\simple_mod_effs_coding_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich\\simple_mod_effs_coding_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich\\simple_mod_effs_coding_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich\\simple_mod_effs_coding_3.csv', comment.char='#')
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
#          U.1.16.1, U.2.16.1, U.1.17.1, U.2.17.1, U.1.18.1, U.2.18.1,
#          #trt_type linear slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.2, U.2.1.2, U.1.2.2, U.2.2.2, U.1.3.2, U.2.3.2, U.1.4.2, U.2.4.2, U.1.5.2, U.2.5.2,
#          U.1.6.2, U.2.6.2, U.1.7.2, U.2.7.2, U.1.8.2, U.2.8.2, U.1.9.2, U.2.9.2, U.1.10.2, U.2.10.2,
#          U.1.11.2, U.2.11.2, U.1.12.2, U.2.12.2, U.1.13.2, U.2.13.2, U.1.14.2, U.2.14.2, U.1.15.2, U.2.15.2,
#          U.1.16.2, U.2.16.2, U.1.17.2, U.2.17.2, U.1.18.2, U.2.18.2,
#          #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.3, U.2.1.3, U.1.2.3, U.2.2.3, U.1.3.3, U.2.3.3, U.1.4.3, U.2.4.3, U.1.5.3, U.2.5.3,
#          U.1.6.3, U.2.6.3, U.1.7.3, U.2.7.3, U.1.8.3, U.2.8.3, U.1.9.3, U.2.9.3, U.1.10.3, U.2.10.3,
#          U.1.11.3, U.2.11.3, U.1.12.3, U.2.12.3, U.1.13.3, U.2.13.3, U.1.14.3, U.2.14.3, U.1.15.3, U.2.15.3,
#          U.1.16.3, U.2.16.3, U.1.17.3, U.2.17.3, U.1.18.3, U.2.18.3,
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
# write.csv(chainsCommunity2, 'bayesian_output_summary_nointeractions_12182017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_nointeractions_12182017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,3254:5725]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 3254:5725])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,3254:5725]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 3254:5725])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_nointeractions_12182017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_nointeractions_12182017.csv')

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
  # mutate(color=ifelse(trt_type=='CO2', '#fabebe44', ifelse(trt_type=='N', '#3cb44b44', ifelse(trt_type=='P', '#aa6e2844', ifelse(trt_type=='drought', '#ffe11944', ifelse(trt_type=='irr', '#0082c844', ifelse(trt_type=='precip_vari', '#46f0f044', ifelse(trt_type=='burn', '#f5823144', ifelse(trt_type=='mow_clip', '#00808044', ifelse(trt_type=='herb_rem', '#f032e644', ifelse(trt_type=='temp', '#e6194b44', ifelse(trt_type=='plant_mani', '#911e6444', ifelse(trt_type=='R*burn', '#f5823144', ifelse(trt_type=='R*mow_clip', '#00808044', ifelse(trt_type=='R*herb_rem', '#f032e644', ifelse(trt_type=='R*temp', '#e6194b44', ifelse(trt_type=='R*plant_mani', '#911e6444', ifelse(trt_type=='R*R', '#80800044', ifelse(trt_type=='N*N', '#80000044', ifelse(trt_type=='all_resource', '#00008044', '#80808044'))))))))))))))))))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1721255)+(0.3228545)}, size=2, xlim=c(0,',
                       '*x^2)*(0.2468055)+(-0.07086403)}, size=2, xlim=c(0,'),
         curve5='), colour=grey) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_nointeractions_12182017.csv', row.names=F)

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
  stat_function(fun=function(x){(-0.44814002034765 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6737189045605 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.57656166673475 + 0.2373468476496*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9535768582321 + 0.2800205370192*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.16230827855 + 0.29697051301605*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.58428166790595 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.451318978352735 + 0.222370893288315*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2913263159398*x + -0.0157383046024025*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2605746633379*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.84816148797 + 0.19006941349798*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.75152184955 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.10134378235 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.746146405235 + 0.17809386812904*x + -0.0100788991943695*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.80563461916425 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.65169187463275 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.59149852130795 + 0.17622684328065*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.744970629117 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.49534005506134 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.60100709888945 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6742947203706 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7082033363745 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.698931243132 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.28577200495 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7511662987145 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.33906652435 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.2130292084 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.2642490971 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.60071673629125 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.94319199024 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8302979922 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7707234865165 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9631452608 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.270101737223245*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7623730804 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6865080153965 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.58975419952605 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.808328366975 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.00642144355 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.88607152908 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.2314326165 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.1230741676 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.2897475861 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.12687686355 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.09330397075 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7099608623525 + 0.20276517783835*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.648708078018 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.556323738185 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4508155697646 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.791544743805 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.682370587442 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7814392574505 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.0163801622 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.04337934655 + 0.1443832485197*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.82846354776 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.94227405656 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.824122753025 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.76040537237 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.56734558777755 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.86827261392 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.77396553641 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6029858483745 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.737228869646 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1727292667631*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.20888244671405*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.159275446271075*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.19061995636935*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.697375521325 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.734180549523 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8394972646057 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.01111492757 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.0626435601 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.11923801365 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.84427914446 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.577208691996285 + 0.266099879540755*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.682810433457 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.24199465978363*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5491467958168 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.55705235078489 + 0.3536221611082*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.46800184594895 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.675619878027 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7607116057005 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.10723526345 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.11998570185 + 0.176406807473235*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.93777878265 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.10148345605 + 0.16370313893315*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.04877881315 + 0.170773663284285*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.2832893885 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.06040474965 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.06239390265 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.926898237335 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.446093667494 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-1.3449061373 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.54344250681025 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6461666863153 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.623075087778 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.701097557016 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.35872648 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.61388763083325 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.91815863525 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.66735217255 + 0.20414414245339*x + -0.0109114004430295*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.46063113140395 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9813246278 + 0.202384792904303*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.90536963977 + 0.21235616819475*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.938143013501 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.69505114343015 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7553887020155 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0301154544 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.580641400992 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5554804950385 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.67899540821435 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7799299371295 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.473282507053775 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.677657906297 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.443140619639585 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4565902103268 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.619460584068 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.17370706076561*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8095933252645 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.547590794003745 + 0.133690175658795*x + -0.01042325004885*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8706287579 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.913893087285 + 0.201355772824904*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.3960826697 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.034771024 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6067123567765 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.131037286725297*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.23129879762315*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.872464457595 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.96072514245 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1474269091928*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2492604688805*x + -0.013094425191*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5995455899195 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7992203474 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.92394462215 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.17904424155 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.948356778355 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.806620631377 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9143859799 + 0.3492171241555*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.86482287356 + 0.23087770509755*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.610222683072 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.823947914585239 + 0.17656123208276*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.77693952215775 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.467304448559306 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.45331849422005 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.78310080897 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2592088758325*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4161527788039 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6041021085019 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.492843481222965 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.532725366751045 + 0.42926502337685*x + -0.061008418579975*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.980673236946 + 0.4889566964432*x + -0.06306428805211*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5099323700073 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.534171646398 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5227292499093 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.23326295512113*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6624088242565 + 0.2641450253441*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.60771994673305 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.04202428465 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.1928420583 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6754794509171 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.224775850830435*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.264261967510735*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6483523577055 + 0.27068437991*x + -0.0161814137554*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.881085052355 + 0.15158982638195*x + -0.0110539824308965*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6688324423377 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2542779047371*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.733813522865 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6213912818327 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.770048759749 + 0.25785036798205*x + -0.013883975022112*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.29582700845 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.395474504 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.21837487585 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.60793305143535 + 0.365749719805*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.55818032145031 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.57601978825945 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7377668132455 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6525815214157 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.61777666108805 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.720678431832 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.854324467513 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.98068191065 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.1317383665 + 0.241202433529125*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.14377021155 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.2152648928 + 0.3519971709215*x + -0.016533557753766*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.063288847845 + 0.271440072758285*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.611426267186 + 0.2416057379655*x + -0.0151247910471315*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7924539746536 + 0.269300070799*x + -0.01658935226487*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.64617258053335 + 0.21027808455475*x + -0.011555171735874*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5343058579813 + 0.14174945581575*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.67304604921335 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.331451778776825*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5281556765121 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.17640318881005*x + -0.0067102525467835*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2408796312795*x + -0.01113524060675*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.600644536116375 + 0.153647839120555*x + -0.00894299217134*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.5543865236099 + 0.2442711781369*x + -0.014562123444*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.154134316346255*x + -0.00737748033751355*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.143422126511785*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.513156647914955 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.0248222361 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.21942592425105*x + -0.016259238507875*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.23476494845611*x + -0.0138645090435055*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.23295509736695*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8766263366995 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.29685923591*x + -0.012807072556327*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2839786325068*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.72869448765185 + 0.34269711135*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.28002857265 + 0.40965565518*x + -0.0153406853815139*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.221616859 + 0.4031837058*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.0769545476615 + 0.30998139546413*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.17598085165 + 0.3371834777261*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.6570347739873 + 0.25364421466001*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.53937648744915 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.467280689043 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.55345976166245 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.786873039904 + 0.403837181215*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0547270007 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0826224069 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.03755403525 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.962313511985 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.8651494426739 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.87525225504 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.51003912361125 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.51187967036115 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.5466619615181 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.3203586864 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.09893744635 + 0.24642567408175*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.90159496492 + 0.15348108146016*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6633972222104 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.317185672605*x + -0.014552525627135*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.47216982859599 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.214895876512186*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.50650812160799 + 0.18640322030085*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.51216599422535 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.739872104506 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.65600238270335 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.40389569135 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.47266740335 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.5303123985 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.58564406305 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.5564913866 + 0.16913249262254*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.39460035135 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6184288480611 + 0.138705873395*x + -0.00223991987419457*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8316640056 + 0.20504575471*x + -0.0046526207536*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8462797684 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.33113320065*x + -0.0079439597265*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.40983615045*x + -0.009843566309*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.45385532034999 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.14133778365 + 0.089199200906775*x + -0.003311599631235*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.9235542271 + 0.1059356786035*x + -0.0022751933340485*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.00319934745285422*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.7392595128165 + 0.091593599542725*x + -0.00495752627335*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.698017588775 + 0.236145048045*x + -0.0100408590285*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6607201168376 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.50300407864964 + 0.263420737736*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.28194508466*x + -0.013491879955815*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.55721515735*x + -0.0226033336247*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9143653764511 + 0.553405985*x + -0.02688931900565*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.71940268434 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.354839792965*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5091867731*x + -0.020042494334853*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.90720216809 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.978223981637 + 0.24248983555029*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6262381792182 + 0.315007243505*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.42248339679974 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.420790440117565 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.45035770447075 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.42065317054714 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.48860914535775 + 0.156811158489465*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2151747318343*x + -0.015227241684025*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.42149756549439 + 0.164627647537461*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.415241317778085 + 0.17933115729117*x + -0.01237806911082*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.24678312857415*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2969163020167*x + -0.0166863550057435*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32988218133*x + -0.014923581381045*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.366095096465*x + -0.02002098934055*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.196122576739505*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.129663524511532*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0066557348359675*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.14782650549635*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1997976383256*x + -0.0078938838838015*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.006462995541729*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.293690716289*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2573156995745*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.5853929758748 + 0.24419268870739*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.83679992683 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.1520808433 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36206040138*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.698808709357 + 0.3275918492125*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7121287456395 + 0.394646640695*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.00885885385 + 0.363072423345*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.275905560004*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2932583030805*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.561810771655*x + -0.0625193182934801*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.563317390695*x + -0.06236207255555*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.862723152122 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7060741268715 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.468534709594005 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.212560494205665*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.70405840146975 + 0.522905063163*x + -0.060898558411291*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.48944390204*x + -0.059530798213179*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.45471947203533 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1721255)+(0.3228545)}, size=1.5, xlim=c(0,2), colour='grey') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.485664000 + 0.154460000*x + -0.009882655*x^2)*(0.1721255)+(0.3228545)}, size=5, xlim=c(0,7), colour='black')

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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.562062220920735 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6278256545672 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.35747584558505*x + 0.028146845920683*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.445592150548*x + 0.033983191920355*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.387279730482*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.676771339251 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.77127548616685 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.74135372561865 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.23318849326035*x + 0.0180057330334565*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.21266814042338*x + 0.0200479786349565*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.17006387155 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.654678027753 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.635108664013235 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7197354313638 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7579080068368 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.719484253723 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.67827430580057 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.726142375889605 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.67937238301465 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.47896318959794 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.9237546850805 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.176077994025375*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.688157631615 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.749056621726 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.852613523921 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8435285637245 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.667715770500365 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.8513883501 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.98741500507 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6023554543565*x + 0.04924172067645*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9274184821142 + -0.490364108744*x + 0.040613571768315*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.599228906344*x + 0.0505102650469*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.09706721013 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.287550394659525*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.628430475951195 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.56524701062117 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.70468337328885 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.5207900422259 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.51835870803508 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.666673052788 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6763752478897 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.539993850131452 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.16337098205185*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9923408323 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.63359183388765 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.18530531430191*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.783906888798735 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7194978120985 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.55166637874605 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.5762051627516 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.02993554361 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.91953342609 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.99738030211 + 0.3600703638465*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5405824336041 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.98900317431625 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.753193367024605 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.91342123649 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.15492321854755*x + 0.010290361723085*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.5916297874068 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.75135709978575 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.67443708604185 + 0.242415763002815*x + -0.02677898124049*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7288366647537 + 0.2932278944827*x + -0.02732976085515*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.685550330194225 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.563734846574295 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.15891831821 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.61813824521705 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.47397555147555 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.84318603197 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.54989012863909 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.566544945602926 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.24343902093575*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.62638030622995 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.170556278393942*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.222367813310775*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2066625165418*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.60231574702675 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7446526043875 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.33607447165 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.226741652215 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.2152757165 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.24574713205 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.928835546087 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7717113849311 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.5509507138593 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.709347453874535 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.57668427177825 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5722794822006 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.63330249635335 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.67924943297645 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.567488836182115 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2106268734828*x + 0.0147911420971*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.24104947284385*x + 0.013715624814495*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8491748824935 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.754875977746 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.26999826653675*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.298002063459*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.49520523871473 + -0.2521197807155*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5570339343332 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.50487465421305 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.60190532971865 + 0.36384928998705*x + -0.03665392706275*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.56546170144749 + 0.3015545039043*x + -0.0362872180782*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.58354314625*x + -0.060957067225*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.63144020155*x + -0.065711728685*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9945732209475 + 0.60246682395*x + -0.07344659837*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.5858187072111 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.571337711485215 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.335798091198495*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.04570688875 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.337745301 + 0.64591264316*x + -0.0526707974818155*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.31293357408416*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.805617506906 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5775207995012 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.59051939360345 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.681987566471 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.85399418557 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.56073834294625 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.516299389489465 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.026045195126137*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.60685813514645 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.64352399413385 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.715286278964 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.5333383019 + 0.78360319765*x + -0.077315317805*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6732318688665*x + -0.06924341108*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3357087143273*x + -0.04833686399305*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.51703562241065 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.07537483753392*x + 0.00274937431974765*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.2527741039 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.7332147123 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.91412218899 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.510264279540055 + -0.30693711075*x + 0.0079323530695*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.23484257534*x + 0.0056719435766*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.916481461585 + -0.1395081436335*x + 0.0042799653413765*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-1.0967199721 + -0.12642922889145*x + 0.0033077773005295*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1778068312325*x + 0.00475609631154865*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2757301786692*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.842031442*x + 0.0671733462*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.9580685211*x + 0.078705406165*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8147746432711 + -0.4060826734345*x + 0.0389658590254*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.73757163815*x + 0.059119805255*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.8485492098*x + 0.061967349395*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.41517478197*x + 0.030663789116455*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.58105932485*x + 0.041580783608*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.5348173200833 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.68995539114978 + -0.24306153802976*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.383740460635*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.48005893044174 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.5071745354483 + 0*x + -0.014114747624125*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.530296782585565 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.167653442444815*x + 0.01069458640744*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.20424086329665*x + 0.0106970914524975*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.20379826618521*x + 0.013586091904945*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.42966136492*x + 0.02675542729*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.551755319501395 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.63917486916065 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5379542895225*x + -0.0462506546532*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(1.79842718915 + 0*x + -0.0242185806422*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.22107148012 + 0*x + -0.023960086469533*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.2870104584 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.560483686432125 + 0*x + -0.017151481919433*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.248871096073565*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5095164581934 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.92944928648 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.63088607469 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2468055)+(-0.07086403)}, size=1.5, xlim=c(0,2), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.259751500 + 0*x + 0*x^2)*0.2468055 + -0.07086403}, size=5, xlim=c(0,7), colour='black')

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
#          U.1.16.1, U.2.16.1, U.1.17.1, U.2.17.1, U.1.18.1, U.2.18.1,
#          #trt_type linear slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.2, U.2.1.2, U.1.2.2, U.2.2.2, U.1.3.2, U.2.3.2, U.1.4.2, U.2.4.2, U.1.5.2, U.2.5.2,
#          U.1.6.2, U.2.6.2, U.1.7.2, U.2.7.2, U.1.8.2, U.2.8.2, U.1.9.2, U.2.9.2, U.1.10.2, U.2.10.2,
#          U.1.11.2, U.2.11.2, U.1.12.2, U.2.12.2, U.1.13.2, U.2.13.2, U.1.14.2, U.2.14.2, U.1.15.2, U.2.15.2,
#          U.1.16.2, U.2.16.2, U.1.17.2, U.2.17.2, U.1.18.2, U.2.18.2,
#          #trt_type quadratic slopes: center digit refers to trts and interactions with anpp and gamma diversity
#          U.1.1.3, U.2.1.3, U.1.2.3, U.2.2.3, U.1.3.3, U.2.3.3, U.1.4.3, U.2.4.3, U.1.5.3, U.2.5.3,
#          U.1.6.3, U.2.6.3, U.1.7.3, U.2.7.3, U.1.8.3, U.2.8.3, U.1.9.3, U.2.9.3, U.1.10.3, U.2.10.3,
#          U.1.11.3, U.2.11.3, U.1.12.3, U.2.12.3, U.1.13.3, U.2.13.3, U.1.14.3, U.2.14.3, U.1.15.3, U.2.15.3,
#          U.1.16.3, U.2.16.3, U.1.17.3, U.2.17.3, U.1.18.3, U.2.18.3,
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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_nointeraction_12182017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_nointeraction_12182017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  mutate(type=paste(predictor2, parameter, sep='_'))



###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor2!='trt_type'), aes(x=type, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.4)) +
  scale_y_continuous(limits=c(-0.8, 0.4), breaks=seq(-0.5, 0.5, 0.5)) +
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
  scale_y_continuous(limits=c(-0.8, 0.4), breaks=seq(-0.5, 0.5, 0.5)) +
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


#NOTE: START HERE - need to make the treatment figures for mean and richness (replacing plotmani, anpp, gamma, map, mat)
###treatment, experiment, and site level driver effects from bayesian output (Figure 2)
#mean plots --------------------------------------------------------
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='mean'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2', 'overall'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=8.5), linetype='dashed') +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='mean'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2', 'overall'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=8.5), linetype='dashed') +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='mean'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2', 'overall'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=8.5), linetype='dashed') +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#richness plots --------------------------------------------------------
richnessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='richness'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2', 'overall'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=8.5), linetype='dashed') +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

richnessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='richness'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2', 'overall'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=8.5), linetype='dashed') +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

richnessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='richness'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2', 'overall'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations', 'Overall')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=8.5), linetype='dashed') +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.14,0.1)) +
  coord_flip()

#plot all together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,3))) 
print(meanIntPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanSlopePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanQuadPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(richnessIntPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(richnessSlopePlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(richnessQuadPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 3))
#export at 1250x2000




###by magnitude of resource manipulated---------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
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

Ndispersion <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\n_dispersion_posteriors.csv', comment.char='#')
NdispersionMean <- as.data.frame(colMeans(Ndispersion))%>%
  add_rownames('parameter')
names(NdispersionMean)[names(NdispersionMean) == 'colMeans(Ndispersion)'] <- 'mean'
NdispersionSD <- as.data.frame(colSd(Ndispersion))%>%
  add_rownames('parameter')
names(NdispersionSD)[names(NdispersionSD) == 'colSd(Ndispersion)'] <- 'sd'
NdispersionOverall <- NdispersionMean%>%
  left_join(NdispersionSD)

dispersionNPlotFinal <- ggplot(data=subset(rawTrt, n>0&plot_mani==1), aes(x=n, y=dispersion_change)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(name='Dispersion Difference') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.4, label='(d)', size=12, hjust='left')


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

Nevenness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\n_simpson_posteriors.csv', comment.char='#')
NevennessMean <- as.data.frame(colMeans(Nevenness))%>%
  add_rownames('parameter')
names(NevennessMean)[names(NevennessMean) == 'colMeans(Nevenness)'] <- 'mean'
NevennessSD <- as.data.frame(colSd(Nevenness))%>%
  add_rownames('parameter')
names(NevennessSD)[names(NevennessSD) == 'colSd(Nevenness)'] <- 'sd'
NevennessOverall <- NevennessMean%>%
  left_join(NevennessSD)

evennessNPlotFinal <- ggplot(data=subset(rawTrt, n>0&plot_mani==1), aes(x=n, y=SimpEven_change)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(name='Evenness Difference') +
  xlab('') +
  annotate('text', x=0.4, y=0.6, label='(b)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(evennessNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600

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
  scale_y_continuous(name='Overall Community Difference') +
  stat_function(fun=function(x){(0.1820251 + 0.0002544999*x)}, size=5) +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.65, label='(c)', size=12, hjust='left')

H2Odispersion <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\h20_dispersion_posteriors.csv', comment.char='#')
H2OdispersionMean <- as.data.frame(colMeans(H2Odispersion))%>%
  add_rownames('parameter')
names(H2OdispersionMean)[names(H2OdispersionMean) == 'colMeans(H2Odispersion)'] <- 'mean'
H2OdispersionSD <- as.data.frame(colSd(H2Odispersion))%>%
  add_rownames('parameter')
names(H2OdispersionSD)[names(H2OdispersionSD) == 'colSd(H2Odispersion)'] <- 'sd'
H2OdispersionOverall <- H2OdispersionMean%>%
  left_join(H2OdispersionSD)

dispersionH2OPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=dispersion_change)) +
  geom_point(size=5) +
  scale_y_continuous(name='Dispersion Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.3, label='(d)', size=12, hjust='left')

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
  scale_y_continuous(name='Richness Difference') +
  xlab('') +
  annotate('text', x=-80, y=0.6, label='(a)', size=12, hjust='left')

H2Oevenness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\h20_simpson_posteriors.csv', comment.char='#')
H2OevennessMean <- as.data.frame(colMeans(H2Oevenness))%>%
  add_rownames('parameter')
names(H2OevennessMean)[names(H2OevennessMean) == 'colMeans(H2Oevenness)'] <- 'mean'
H2OevennessSD <- as.data.frame(colSd(H2Oevenness))%>%
  add_rownames('parameter')
names(H2OevennessSD)[names(H2OevennessSD) == 'colSd(H2Oevenness)'] <- 'sd'
H2OevennessOverall <- H2OevennessMean%>%
  left_join(H2OevennessSD)

evennessH2OPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=SimpEven_change)) +
  geom_point(size=5) +
  scale_y_continuous(name='Evenness Difference') +
  xlab('') +
  annotate('text', x=-80, y=0.3, label='(b)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(dispersionH2OPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(evennessH2OPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessH2OPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(meanH2OPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600