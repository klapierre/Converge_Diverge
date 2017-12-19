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

#only run to generate initial chains files
#raw chains data --------------------------------------------------------
memory.limit(size=50000)
chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_ints\\simple_mod_effs_coding_interactions_3.csv', comment.char='#')
chains4 <- chains4[-1:-5000,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)


#density plot of chains --------------------------------------------------------
plot(density(chainsCommunity$effs_G.1.2.68))
plot(density(chainsCommunity$mu.1.2))
plot(density(chainsCommunity$mu.1.3))


#get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
#mean change are the 1's, richness are the 2's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__,
         #trt_type intercepts (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1.1.1, U.2.1.1,
         U.1.2.1, U.2.2.1,
         U.1.3.1, U.2.3.1,
         U.1.4.1, U.2.4.1,
         U.1.5.1, U.2.5.1,
         U.1.6.1, U.2.6.1,
         U.1.7.1, U.2.7.1,
         U.1.8.1, U.2.8.1,
         U.1.9.1, U.2.9.1,
         U.1.10.1, U.2.10.1,
         U.1.11.1, U.2.11.1,
         U.1.12.1, U.2.12.1,
         U.1.13.1, U.2.13.1,
         U.1.14.1, U.2.14.1,
         U.1.15.1, U.2.15.1,
         U.1.16.1, U.2.16.1,
         U.1.17.1, U.2.17.1,
         U.1.18.1, U.2.18.1,
         #trt_type linear slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1.1.2, U.2.1.2,
         U.1.2.2, U.2.2.2,
         U.1.3.2, U.2.3.2,
         U.1.4.2, U.2.4.2,
         U.1.5.2, U.2.5.2,
         U.1.6.2, U.2.6.2,
         U.1.7.2, U.2.7.2,
         U.1.8.2, U.2.8.2,
         U.1.9.2, U.2.9.2,
         U.1.10.2, U.2.10.2,
         U.1.11.2, U.2.11.2,
         U.1.12.2, U.2.12.2,
         U.1.13.2, U.2.13.2,
         U.1.14.2, U.2.14.2,
         U.1.15.2, U.2.15.2,
         U.1.16.2, U.2.16.2,
         U.1.17.2, U.2.17.2,
         U.1.18.2, U.2.18.2,
         #trt_type quad slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1.1.3, U.2.1.3,
         U.1.2.3, U.2.2.3,
         U.1.3.3, U.2.3.3,
         U.1.4.3, U.2.4.3,
         U.1.5.3, U.2.5.3,
         U.1.6.3, U.2.6.3,
         U.1.7.3, U.2.7.3,
         U.1.8.3, U.2.8.3,
         U.1.9.3, U.2.9.3,
         U.1.10.3, U.2.10.3,
         U.1.11.3, U.2.11.3,
         U.1.12.3, U.2.12.3,
         U.1.13.3, U.2.13.3,
         U.1.14.3, U.2.14.3,
         U.1.15.3, U.2.15.3,
         U.1.16.3, U.2.16.3,
         U.1.17.3, U.2.17.3,
         U.1.18.3, U.2.18.3,
         #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
         D.1.1.1, D.2.1.1,
         D.1.1.2, D.2.1.2,
         D.1.1.3, D.2.1.3,
         #richness intercept, linear, and quad slopes (center digit): 2=gamma diversity
         D.1.2.1, D.2.2.1,
         D.1.2.2, D.2.2.2,
         D.1.2.3, D.2.2.3,
         #MAP intercept, linear, and quad slopes (center digit): 1=MAP
         E.1.1.1, E.2.1.1,
         E.1.1.2, E.2.1.2,
         E.1.1.3, E.2.1.3,
         #MAT intercept, linear, and quad slopes (center digit): 2=MAT
         E.1.2.1, E.2.2.1,
         E.1.2.2, E.2.2.2,
         E.1.2.3, E.2.2.3,
         #overall intercept, linear, and quad slopes
         mu.1.1, mu.2.1,
         mu.1.2, mu.2.2,
         mu.1.3, mu.2.3)%>%
  gather(key=parameter, value=value, U.1.1.1:mu.2.3)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))

write.csv(chainsCommunity2, 'bayesian_output_summary_final plots_11292017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_final plots_11292017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,6002:8449]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 6002:8449])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,6002:8449]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 6002:8449])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_11292017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_11292017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=8, 7, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1405501)+(0.2906267), (intercept+linear*7+quadratic*7^2)*(0.2254782)+(-0.03185368)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1405501)+(0.2906267),
                    (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2254782)+(-0.03185368)))%>%
  mutate(color=ifelse(trt_type=='CO2', '#fabebe44', ifelse(trt_type=='N', '#3cb44b44', ifelse(trt_type=='P', '#aa6e2844', ifelse(trt_type=='drought', '#ffe11944', ifelse(trt_type=='irr', '#0082c844', ifelse(trt_type=='precip_vari', '#46f0f044', ifelse(trt_type=='burn', '#f5823144', ifelse(trt_type=='mow_clip', '#00808044', ifelse(trt_type=='herb_rem', '#f032e644', ifelse(trt_type=='temp', '#e6194b44', ifelse(trt_type=='plant_mani', '#911e6444', ifelse(trt_type=='R*burn', '#f5823144', ifelse(trt_type=='R*mow_clip', '#00808044', ifelse(trt_type=='R*herb_rem', '#f032e644', ifelse(trt_type=='R*temp', '#e6194b44', ifelse(trt_type=='R*plant_mani', '#911e6444', ifelse(trt_type=='R*R', '#80800044', ifelse(trt_type=='N*N', '#80000044', ifelse(trt_type=='all_resource', '#00008044', '#80808044'))))))))))))))))))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,',
                       '*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,'),
         curve5='), colour=',
         curve6=') +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_11292017.csv', row.names=F)


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
  annotate('text', x=0, y=1, label='(c)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.4012628709935*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.7583679733232 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0.8071464847643 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(1.25632954507 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(1.3368254263 + 0.44065249811175*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0.44891016030665*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0.76412088371765 + 0.433283498395045*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.4123390246769*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(-0.9187456813235 + 0.35771607722442*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(-0.76460730967295 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.2198929882 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.65388850368835 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.62769873574235 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.62433595165715 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(-0.6541027631952 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.858137725339 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.83330402161235 + 0.48026275403255*x + -0.057680249833725*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#911e6444') +
  stat_function(fun=function(x){(-1.3572157492 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#911e6444') +
  stat_function(fun=function(x){(-0.8676565797835 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#911e6444') +
  stat_function(fun=function(x){(-1.44608935405 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.2670738015 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.4400653977 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0.432995144111975*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0.344716917315*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(-0.748981960765 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.08257416431 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-0.82272000595175 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(-0.89031212342495 + 0.48859850885995*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-1.09299461115 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0.66346400928475*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-0.81555105942915 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-0.7054134957541 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-0.6102024636566 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-0.87697174199925 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(-0.9757471753935 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#46f0f044') +
  stat_function(fun=function(x){(-0.867333997734 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#46f0f044') +
  stat_function(fun=function(x){(-1.2087106213 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#46f0f044') +
  stat_function(fun=function(x){(-1.07718817821 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#ffe11944') +
  stat_function(fun=function(x){(-1.298873575975 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.2336969325 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.188173019 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.7804997118885 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(-0.54850636801455 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(-0.8177876932445 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.6625480932235 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(-0.59887714879935 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(-1.09013964713 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(-0.879488513495 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.8782318870955 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(-0.7919549730985 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(-0.8990846918 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(-0.642028982793945 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.7830868722655 + 0.48138756007515*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(-0.717611359869 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(-0.70905036201935 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(-0.50436135082859 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.367506885716935*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + 0.32195207573494*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.6739305351784 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.694581461771663 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.812494574388 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.9960699663535 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-1.03590653994 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-1.16374085215 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0.37944385680615*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.998480821655 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.544799406647365 + 0.352287280108845*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.62315867327795 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.8062823147775 + 0.5698880334855*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#fabebe44') +
  stat_function(fun=function(x){(-0.5579581690532 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.836755308601 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.07914435392 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.18871094975 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80808044') +
  stat_function(fun=function(x){(-0.90373968073952 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(-1.050638414345 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.11691503805 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(-1.22371821005 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.07988138718 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(-0.919343578711135 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.913233479287 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.39125073485 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80800044') +
  stat_function(fun=function(x){(-0.6606629109098 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.554741876914935 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.6633467486559 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.379633202896695*x + -0.0758514285117285*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.77802188753 + 0.315780601444965*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.784362904928 + 0.5049757240435*x + -0.0748276534818*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-1.46289440055 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.085898166 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-1.072291331485 + 0.3213925129511*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.958930267735 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.79232282827575 + 0.41121955167485*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.7220686307082 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.5200699873755*x + -0.06116813123296*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.046472193645 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.611133370277485 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.54198320701935 + 0.390716069731955*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.339934309399805*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.8313636327555 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.90623676721 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.66559029816065 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.61975339760965*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(-0.8551169766165 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0.406979127293415*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(-0.5569435532165 + 0.3540776671323*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(-0.7388005126781 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(-0.8145522675815 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(-1.3658701653 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(-1.09688888865 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.6269728219431 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.42511716281295*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-0.8172112311105 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-0.774965226291 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.6410966222815 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.777660680527 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(-0.796823748868 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.082213522425 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(-0.7660428389134 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.63489347589775 + 0*x + 0.0523220225153485*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.61859947266595 + 0*x + 0.04973158048488*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(-0.60710094529415 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(-0.7929607228801 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f032e644') +
  stat_function(fun=function(x){(-0.708408540896075 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0.52061364446275*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(-0.78256855621315 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.5869547525342*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(1.2943589015 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(1.1204098234 + 0.5593469662435*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(1.6901135689 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0.8430069644665 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0.625053812936005 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0.54079542222155 + 0.49450630304955*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.470321183020605*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0.5626734663698 + 0.6577424840999*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0.406854882725565*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0.48778183003775*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(-0.871994439702 + 0.46702430489135*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(-0.889166195202 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.7221835486451 + 0.4232102013892*x + -0.05816032376265*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#46f0f044') +
  stat_function(fun=function(x){(-0.583973126531045 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.784347553924 + 0.463970133060035*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.2747057733 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(-1.39646886105 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.35040394435 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-0.599609939319 + 0.3975753132252*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.51136256980555 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-0.5151780234359 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.72943426405265 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.61683347491355 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-0.71760802135335 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-0.66587656256195 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.9218173321655 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(-1.19378430035 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(-1.052801199265 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.08087409217 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(-0.5891042850593 + 0.31993374188177*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.6264978771912 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80808044') +
  stat_function(fun=function(x){(-0.64130630289694 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(-0.7710726548685 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0.398118168346195*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0.63612081802245 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.74363002350595 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(-0.71634044144855 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(-0.6470840041293 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(-0.968342163935 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0.5918792131311*x + -0.0590731801206605*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.4790009846625*x + -0.055341473720835*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0.6006746152125*x + -0.0505930197448195*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80000044') +
  stat_function(fun=function(x){(-0.643021549393 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.409164765236775*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.505471393931*x + -0.048527998777715*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.46564203131675*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0.40057499656915*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(-1.1171813443 + 0*x + 0.0497717671414875*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(-1.1632851733 + 0*x + 0.04740538931315*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.163115113955 + 0.34765628566555*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.8402072519925 + 0*x + 0.0608086760739*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#0082c844') +
  stat_function(fun=function(x){(0.71464157065925 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(-0.55621331946745 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0.7627616714025 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.777227509482 + 0.433081790366*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(-1.09012004183 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(-1.01737812342 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.95671731399 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(-0.8814877267085 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#aa6e2844') +
  stat_function(fun=function(x){(-0.80032028543595 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.6661806230727 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(-1.27657391165 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#f032e644') +
  stat_function(fun=function(x){(-0.968819193355 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(-0.8031922361332 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.612844459387415 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0.32036166353851*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.66311511556173 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.6438654389834 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(-1.633670355 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.71078433525 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#0082c844') +
  stat_function(fun=function(x){(-1.6589288663 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#0082c844') +
  stat_function(fun=function(x){(-1.7352721836 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.64755338185 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(-1.7198473641 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.71831910562*x + -0.06589346494375*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0.8056754049766 + 0.66510733798*x + -0.06587146584845*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-1.07809485702 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.7769019068612 + 0.388385247677555*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(-0.6483154332797 + 0.57731159769*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(-0.57700710845705 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0.9430334101*x + -0.0947155232775*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.9795988429*x + -0.0887455360265*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0.5246403620664*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0.85397946105*x + -0.050119225409545*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#00808044') +
  stat_function(fun=function(x){(-0.574247286603661 + 0.4868187294485*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.42168547703591*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.8408322382*x + -0.04330722281914*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(0.7921824355055 + 1.16228727945*x + -0.09115148191775*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(-0.5635033800755 + 0.566279898615*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0.81465047705*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(-1.130349931865 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(-1.04665390148 + 0.351327234423723*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(-0.918157971821 + 0.6825127068455*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0.32430596972495*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.3190661823243*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0.38685524317445*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.3997235346712*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0.5467450541785*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0.5985102724485*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0.3774931497332*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#80000044') +
  stat_function(fun=function(x){(-0.5186719819985 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(-0.648178089112945 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#80808044') +
  stat_function(fun=function(x){(-0.757606268484 + 0.44045118730675*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(-1.14582915447 + 0.3870469901497*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.5648745273085*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(-0.762881827972 + 0.36604007698886*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(-1.13538856001 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0.33193527661562*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0.5733516837505*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0.55804684366165*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.536529141593195*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-1.077438885505 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(-0.618672106862215 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0.8060571794465 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0.41862525441165*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.48749166199625*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.41529836500925*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1405501)+(0.2906267)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.46463900 + 0.16878700*x + -0.01217110*x^2)*(0.1417218)+(0.2924074)}, size=5, xlim=c(0,7), colour='black')

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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(-0.622854245701425 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0.73405918526705 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.72311594470725 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#911e6444') +
  stat_function(fun=function(x){(0.653259704558765 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(1.04458147801 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0.66503953391175 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + -0.743347892943831*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0.68381256085053 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0.874767096193 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0.45523254811895*x + -0.07700026552051*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0.7119606146838 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + 0*x + -0.064180445847425*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.86408324046545 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.8968912574735 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.68321236756195 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(1.71738185635 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(1.2508403282925 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + -0.921566612057*x + 0.0896461540797755*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.8963797535084 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + -0.78568744757065*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0.6719371457996 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + -0.4994213635388*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + -0.44011625852995*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.66813643073235*x + -0.090223716802477*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.7168133306936 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.6259315423296*x + -0.08036126545458*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-1.0389476338265 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0.60020794191615*x + -0.07917460174415*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.642975067852349 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-1.70973999325 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-0.925306545239 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-1.22266711174 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(-1.237639992427 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0.7860655466401 + 0.66587232630795*x + -0.077225091241855*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0.6915459047842 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0.7252683087843*x + -0.09118949352175*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 1.0711170101*x + -0.144224576645*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 1.03152963845*x + -0.13591999606*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0.7255810233155*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0.779024239049 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(-0.6756184682076 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + -0.472030685445*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.675672102293375 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + -0.47979243934395*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0.641327223173425 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(1.3955583152 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(1.10918875744 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(1.298352816725 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(1.22725115866 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0.907297949348 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.8198059793149 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + -0.73313945030475*x + 0.08528495709356*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0.643296804858746*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#46f0f044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0.7133798610328*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.7643087392097*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.61473767531615*x + -0.0874410263211005*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80000044') +
  stat_function(fun=function(x){(1.02864691439215 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#f5823144') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.5895219970617*x + -0.10879135402425*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0.58240083851905*x + -0.120135815405*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 1.0716925986*x + -0.147246742935*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 1.0496329973*x + -0.14097552918*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0.99989801095*x + -0.13179481546*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0.744622812377*x + -0.1073207887645*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#00008044') +
  stat_function(fun=function(x){(0 + 1.2153348659*x + -0.17579276086*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#0082c844') +
  stat_function(fun=function(x){(0.762319256562645 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(1.43001568575 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.92726120909598 + 1.17227245185*x + -0.1457787419635*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0.696300817256204 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#00008044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.445278726402845*x + -0.070936627026012*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#f032e644') +
  stat_function(fun=function(x){(0.728307352795 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0.9969422487999 + 1.6559920136*x + -0.20438948898*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0.8082855327953*x + -0.09986151122535*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#80800044') +
  stat_function(fun=function(x){(0.80631157527635 + -0.648762881505*x + 0.0632989725450055*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + -1.4583152831*x + 0.143078047928*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + -1.56381473605*x + 0.15281229253*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0.80049704439485 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(1.3580134061 + -1.00431515575*x + 0.083413464454365*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + -0.9708366412*x + 0.0890371565694*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + -0.97461538988*x + 0.09860420724325*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + -1.0113163244*x + 0.10120485637835*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + -0.7432098973375*x + 0.093370419276205*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + -1.0486949018*x + 0.121182950215*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,5), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#00808044') +
  stat_function(fun=function(x){(0 + -1.2305548849*x + 0.1137863887702*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#00808044') +
  stat_function(fun=function(x){(0 + -1.45437781675*x + 0.15353512657*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,3), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + -0.9727205032*x + 0.085512424906992*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + -1.375919728*x + 0.12772744845*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,7), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + -0.7878523232345*x + 0.07518312088252*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#0082c844') +
  stat_function(fun=function(x){(0 + -0.94581513978*x + 0.0964063521375*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + -0.672295544594*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0716171940959515*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + -0.059045912500502*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#ffe11944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0.50186012831725*x + -0.076403165574457*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + -0.062217759354741*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0.60980457421155*x + -0.086994784946565*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + -0.0769828396103195*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + -0.076881751240825*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0.84493910092021*x + -0.13357671787587*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0.475836805345105*x + -0.094887727736725*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#0082c844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + -0.08600064339044*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 1.16599083535*x + -0.1513129008*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#80808044') +
  stat_function(fun=function(x){(0 + 1.17743735115*x + -0.18786606812*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,4), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0.852601612845*x + -0.1335183443645*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0.83075555007055 + 1.385016787*x + -0.21105707395*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 1.06419292735*x + -0.17076401452*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + -0.065200055409795*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80800044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#aa6e2844') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#fabebe44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,6), colour='#e6194b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00008044') +
  stat_function(fun=function(x){(0 + -0.5465666790843*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#911e6444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#00808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#f032e644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0.960675314838 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80000044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#80808044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2254782)+(-0.03185368)}, size=2, xlim=c(0,2), colour='#3cb44b44') +
  #overall line
  stat_function(fun=function(x){(0.29502450 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=5, xlim=c(0,7), colour='black')

# print(richnessPlot) #export at 1200x1000

#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(1,2)))
print(richnessPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
#export at 2400 x 2000


###summary stats from bayesian output --------------------------------------------------------
# # gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(
#              #trt_type intercepts (center digit): trt_type
#              U.1.1.1, U.2.1.1,
#              U.1.2.1, U.2.2.1,
#              U.1.3.1, U.2.3.1,
#              U.1.4.1, U.2.4.1,
#              U.1.5.1, U.2.5.1,
#              U.1.6.1, U.2.6.1,
#              U.1.7.1, U.2.7.1,
#              U.1.8.1, U.2.8.1,
#              U.1.9.1, U.2.9.1,
#              U.1.10.1, U.2.10.1,
#              U.1.11.1, U.2.11.1,
#              U.1.12.1, U.2.12.1,
#              U.1.13.1, U.2.13.1,
#              U.1.14.1, U.2.14.1,
#              U.1.15.1, U.2.15.1,
#              U.1.16.1, U.2.16.1,
#              U.1.17.1, U.2.17.1,
#              U.1.18.1, U.2.18.1,
#              #trt_type linear slopes (center digit)
#              U.1.1.2, U.2.1.2,
#              U.1.2.2, U.2.2.2,
#              U.1.3.2, U.2.3.2,
#              U.1.4.2, U.2.4.2,
#              U.1.5.2, U.2.5.2,
#              U.1.6.2, U.2.6.2,
#              U.1.7.2, U.2.7.2,
#              U.1.8.2, U.2.8.2,
#              U.1.9.2, U.2.9.2,
#              U.1.10.2, U.2.10.2,
#              U.1.11.2, U.2.11.2,
#              U.1.12.2, U.2.12.2,
#              U.1.13.2, U.2.13.2,
#              U.1.14.2, U.2.14.2,
#              U.1.15.2, U.2.15.2,
#              U.1.16.2, U.2.16.2,
#              U.1.17.2, U.2.17.2,
#              U.1.18.2, U.2.18.2,
#              #trt_type quad slopes (center digit)
#              U.1.1.3, U.2.1.3,
#              U.1.2.3, U.2.2.3,
#              U.1.3.3, U.2.3.3,
#              U.1.4.3, U.2.4.3,
#              U.1.5.3, U.2.5.3,
#              U.1.6.3, U.2.6.3,
#              U.1.7.3, U.2.7.3,
#              U.1.8.3, U.2.8.3,
#              U.1.9.3, U.2.9.3,
#              U.1.10.3, U.2.10.3,
#              U.1.11.3, U.2.11.3,
#              U.1.12.3, U.2.12.3,
#              U.1.13.3, U.2.13.3,
#              U.1.14.3, U.2.14.3,
#              U.1.15.3, U.2.15.3,
#              U.1.16.3, U.2.16.3,
#              U.1.17.3, U.2.17.3,
#              U.1.18.3, U.2.18.3,
#              #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
#              D.1.1.1, D.2.1.1,
#              D.1.1.2, D.2.1.2,
#              D.1.1.3, D.2.1.3,
#              #richness intercept, linear, and quad slopes (center digit): 2=gamma diversity
#              D.1.2.1, D.2.2.1,
#              D.1.2.2, D.2.2.2,
#              D.1.2.3, D.2.2.3,
#              #MAP intercept, linear, and quad slopes (center digit): 1=MAP
#              E.1.1.1, E.2.1.1,
#              E.1.1.2, E.2.1.2,
#              E.1.1.3, E.2.1.3,
#              #MAT intercept, linear, and quad slopes (center digit): 2=MAT
#              E.1.2.1, E.2.2.1,
#              E.1.2.2, E.2.2.2,
#              E.1.2.3, E.2.2.3,
#              #overall intercept, linear, and quad slopes
#              mu.1.1, mu.2.1,
#              mu.1.2, mu.2.2,
#              mu.1.3, mu.2.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:mu.2.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(CI=sd*2)%>%
#   separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
#   mutate(parameter=ifelse(level=='mu', predictor, parameter), predictor=ifelse(level=='mu', 'overall', predictor))%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', 'richness'),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          predictor=ifelse(level=='D'&predictor==1, 'ANPP', ifelse(level=='D'&predictor==2, 'rrich', ifelse(level=='E'&predictor==1, 'MAP', ifelse(level=='E'&predictor==2, 'MAT', ifelse(level=='mu', 'overall', 'trt_type'))))))%>%
#   select(level, parameter, variable, predictor, predictor, median, sd, CI)
# 
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_11292017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_11292017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)




###overall responses from bayesian output (Figure S2) --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Mean\nChange') +
  annotate('text', x=3.45, y=-0.8, label='(a)', size=10, hjust='left')

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='dispersion' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.18), breaks=seq(-0.2, 0.2, 0.2)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Dispersion\nChange') +
  annotate('text', x=3.45, y=-0.3, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='richness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.5), breaks=seq(-0.3, 0.5, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Richness\nChange') +
  annotate('text', x=3.45, y=-0.3, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.4, 0.25), breaks=seq(-0.3, 0.3, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.4, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 2400x500




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

#dispersion plots --------------------------------------------------------
dispersionIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='dispersion'), aes(x=predictor, y=median)) +
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

dispersionSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='dispersion'), aes(x=predictor, y=median)) +
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

dispersionQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='dispersion'), aes(x=predictor, y=median)) +
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

#evenness plots --------------------------------------------------------
evennessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='evenness'), aes(x=predictor, y=median)) +
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

evennessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='evenness'), aes(x=predictor, y=median)) +
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

evennessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='evenness'), aes(x=predictor, y=median)) +
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

#plot all together --------------------------------------------------------
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
#export at 2500x2000




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






###comparing different resource manipulation types (Figure 3)
#get means and sd to backtransform
trtInteractionsData <- read.csv('treatment interactions_11152017.csv')%>%
  summarize(mean_mean_change=mean(mean_change), sd_mean_change=sd(mean_change), mean_dispersion_change=mean(dispersion_change), sd_dispersion_change=sd(dispersion_change), mean_richness_change=mean(S_PC), sd_richness_change=sd(S_PC), mean_evenness_change=mean(SimpEven_change), sd_evenness_change=sd(SimpEven_change))

#get counts
trtInteractionsCount <- read.csv('treatment interactions_11152017.csv')%>%
  group_by(trt_type)%>%
  summarise(count=length(mean_change))%>%
  ungroup()


#mean change
trtmean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_mean_posteriors.csv', comment.char='#')*0.1645314+0.3341639
trtmeanMean <- as.data.frame(colMeans(trtmean))%>%
  add_rownames('parameter')
names(trtmeanMean)[names(trtmeanMean) == 'colMeans(trtmean)'] <- 'mean'
trtmeanSD <- as.data.frame(colSd(trtmean))%>%
  add_rownames('parameter')
names(trtmeanSD)[names(trtmeanSD) == 'colSd(trtmean)'] <- 'sd'
trtmeanOverall <- trtmeanMean%>%
  left_join(trtmeanSD)

meantrtPlotFinal <- ggplot(data=subset(trtmeanOverall, parameter!='CO2...Irr'&parameter!='P...K'&parameter!='N...CO2...Irr'), aes(x=parameter, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_y_continuous(name='Overall Community Difference') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr', 'Other', 'Res...Other'), labels=c(expression(paste(CO[2], '(8)')), 'Drought (13)', 'Irrigation (28)', 'N (72)', 'P (22)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (33)', 'N*P*K (18)', 'N*P*K*Irr (2)', 'other (81)', 'R*other (133)')) +
  annotate('text', x=0.5, y=0.6, label='(c)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#dispersion
trtdispersion <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_disp_posteriors.csv', comment.char='#')*0.09661412+0.004211247
trtdispersionMean <- as.data.frame(colMeans(trtdispersion))%>%
  add_rownames('parameter')
names(trtdispersionMean)[names(trtdispersionMean) == 'colMeans(trtdispersion)'] <- 'mean'
trtdispersionSD <- as.data.frame(colSd(trtdispersion))%>%
  add_rownames('parameter')
names(trtdispersionSD)[names(trtdispersionSD) == 'colSd(trtdispersion)'] <- 'sd'
trtdispersionOverall <- trtdispersionMean%>%
  left_join(trtdispersionSD)

dispersiontrtPlotFinal <- ggplot(data=subset(trtdispersionOverall, parameter!='CO2...Irr'&parameter!='P...K'&parameter!='N...CO2...Irr'), aes(x=parameter, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_y_continuous(name='Dispersion Difference') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr', 'Other', 'Res...Other'), labels=c(expression(paste(CO[2], '(8)')), 'Drought (13)', 'Irrigation (28)', 'N (72)', 'P (22)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (33)', 'N*P*K (18)', 'N*P*K*Irr (2)', 'other (81)', 'R*other (133)')) +
  annotate('text', x=0.5, y=0.2, label='(d)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#richness
trtrichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_richness_posteriors.csv', comment.char='#')*0.2518764-0.06504351
trtrichnessMean <- as.data.frame(colMeans(trtrichness))%>%
  add_rownames('parameter')
names(trtrichnessMean)[names(trtrichnessMean) == 'colMeans(trtrichness)'] <- 'mean'
trtrichnessSD <- as.data.frame(colSd(trtrichness))%>%
  add_rownames('parameter')
names(trtrichnessSD)[names(trtrichnessSD) == 'colSd(trtrichness)'] <- 'sd'
trtrichnessOverall <- trtrichnessMean%>%
  left_join(trtrichnessSD)

richnesstrtPlotFinal <- ggplot(data=subset(trtrichnessOverall, parameter!='CO2...Irr'&parameter!='P...K'&parameter!='N...CO2...Irr'), aes(x=parameter, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_y_continuous(name='Richness Difference') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr', 'Other', 'Res...Other'), labels=c(expression(paste(CO[2], '(8)')), 'Drought (13)', 'Irrigation (28)', 'N (72)', 'P (22)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (33)', 'N*P*K (18)', 'N*P*K*Irr (2)', 'other (81)', 'R*other (133)')) +
  annotate('text', x=0.5, y=0.5, label='(a)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#evenness
trtevenness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_simp_posteriors.csv', comment.char='#')*0.09968657+0.01221669
trtevennessMean <- as.data.frame(colMeans(trtevenness))%>%
  add_rownames('parameter')
names(trtevennessMean)[names(trtevennessMean) == 'colMeans(trtevenness)'] <- 'mean'
trtevennessSD <- as.data.frame(colSd(trtevenness))%>%
  add_rownames('parameter')
names(trtevennessSD)[names(trtevennessSD) == 'colSd(trtevenness)'] <- 'sd'
trtevennessOverall <- trtevennessMean%>%
  left_join(trtevennessSD)

evennesstrtPlotFinal <- ggplot(data=subset(trtevennessOverall, parameter!='CO2...Irr'&parameter!='P...K'&parameter!='N...CO2...Irr'), aes(x=parameter, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_y_continuous(name='Evenness Difference') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr', 'Other', 'Res...Other'), labels=c(expression(paste(CO[2], '(8)')), 'Drought (13)', 'Irrigation (28)', 'N (72)', 'P (22)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (33)', 'N*P*K (18)', 'N*P*K*Irr (2)', 'other (81)', 'R*other (133)')) +
  annotate('text', x=0.5, y=0.15, label='(b)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

pushViewport(viewport(layout=grid.layout(2,2)))
print(meantrtPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(dispersiontrtPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(evennesstrtPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnesstrtPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 1800 x 2400




###look at five factor manipulations for mean change (Figure S5)------------------------------------------------
#compare any four factor without N to five factor with N - using raw data final year--------------------------
expRawMean <- expRaw%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarise(n=mean(n), herb_removal=mean(herb_removal), plant_mani=mean(plant_mani))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, plot_mani, n, herb_removal, plant_mani)

meanCompare <- subset(rawTrt, project_name=='e001'|project_name=='e002'|site_code=='NIN'|site_code=='TRA')%>%
  left_join(expRawMean, all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))%>%
  mutate(trt=paste(n_mani, plot_mani, sep=':'))%>%
  mutate(keep=ifelse(experiment_length==treatment_year, 1, ifelse(treatment_year==8, 1, 0)))

#plot without N at four factors, with N at five factors--------------------------
ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6&keep==1), variable='mean_change', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.8, 0.2), name='Overall Community Difference') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,0.8)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.37, label='a*', size=10) +
  annotate('text', x=2, y=0.52, label='b*', size=10) +
  annotate('text', x=3, y=0.73, label='c*', size=10)
#export at 900x900