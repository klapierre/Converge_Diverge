library(ggplot2)
library(grid)
library(tidyverse)

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

rawData <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\for resource type analysis_nov2017\\ForAnalysis_allAnalysisLong.csv')

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
# chains1 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_long\\simple_mod_effs_coding_long_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_long\\simple_mod_effs_coding_long_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_long\\simple_mod_effs_coding_long_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\anpp_rich_long\\simple_mod_effs_coding_long_3.csv', comment.char='#')
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
# write.csv(chainsCommunity2, 'bayesian_output_summary_nointeractions_long_12182017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_nointeractions_long_12182017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,1754:2983]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 1754:2983])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,1754:2983]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 1754:2983])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_nointeractions_long_12182017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_nointeractions_long_12182017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=8, 7, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1862561)+(0.3359782), (intercept+linear*7+quadratic*7^2)*(0.2540484)+(-0.0909037)))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1862561)+(0.3359782),
                         (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2540484)+(-0.0909037)))%>%
  # mutate(color=ifelse(trt_type=='CO2', '#fabebe44', ifelse(trt_type=='N', '#3cb44b44', ifelse(trt_type=='P', '#aa6e2844', ifelse(trt_type=='drought', '#ffe11944', ifelse(trt_type=='irr', '#0082c844', ifelse(trt_type=='precip_vari', '#46f0f044', ifelse(trt_type=='burn', '#f5823144', ifelse(trt_type=='mow_clip', '#00808044', ifelse(trt_type=='herb_rem', '#f032e644', ifelse(trt_type=='temp', '#e6194b44', ifelse(trt_type=='plant_mani', '#911e6444', ifelse(trt_type=='R*burn', '#f5823144', ifelse(trt_type=='R*mow_clip', '#00808044', ifelse(trt_type=='R*herb_rem', '#f032e644', ifelse(trt_type=='R*temp', '#e6194b44', ifelse(trt_type=='R*plant_mani', '#911e6444', ifelse(trt_type=='R*R', '#80800044', ifelse(trt_type=='N*N', '#80000044', ifelse(trt_type=='all_resource', '#00008044', '#80808044'))))))))))))))))))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1862561)+(0.3359782)}, size=2, xlim=c(0,',
                       '*x^2)*(0.2540484)+(-0.0909037)}, size=2, xlim=c(0,'),
         curve5='), colour=grey) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, sep=''))%>%
  mutate(trt_overall=ifelse(trt_type=='CO2'|trt_type=='N'|trt_type=='P'|trt_type=='drought'|trt_type=='irr'|trt_type=='precip_vari', 'single_resource', ifelse(trt_type=='burn'|trt_type=='mow_clip'|trt_type=='herb_rem'|trt_type=='temp'|trt_type=='plant_mani', 'single_nonresource', ifelse(trt_type=='all_resource'|trt_type=='both', 'three_way', 'two_way'))))
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_nointeractions_12182017.csv', row.names=F)

###main figure (Figure 1)
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,28), breaks=seq(0,28,2), labels=seq(1,29,2)) +
  ylim(-10,10) +
  xlab('Standardized Experiment Year') +
  ylab('Overall Community Difference') +
  annotate('text', x=0, y=1, label='(b)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines

  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.549363500 + 0.137805500*x + -0.006562025*x^2)*(0.1862561)+(0.3359782)}, size=5, xlim=c(0,28), colour='black')

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
  
  #overall line
  stat_function(fun=function(x){(0.360073500 + 0*x + 0*x^2)*0.2540484 + -0.0909037}, size=5, xlim=c(0,7), colour='black')

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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_nointeraction_long_12182017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_nointeraction_long_12182017.csv')

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