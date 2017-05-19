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

rawData <- read.csv('ForBayesianAnalysis_10yr_Dec2016.csv')

rawData2<- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_disp=mean(dispersion_change), std_disp=sd(dispersion_change), mean_rich=mean(S_PC), std_rich=sd(S_PC), mean_even=mean(SimpEven_change), std_even=sd(SimpEven_change)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani))

expInfo2$study <- seq.int(nrow(expInfo2))
  

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

trtDetail <- expRaw%>%
  mutate(site_code=as.factor(site_code), project_name=as.factor(project_name), community_type=as.factor(community_type), treatment=as.factor(treatment))%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

################################################################################
################################################################################

# #only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)
# 
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsCommunity$mu.1.1))
# plot(density(chainsCommunity$mu.1.2))
# plot(density(chainsCommunity$mu.1.3))
# 
# 
# #get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
# #mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
# chainsCommunity2 <- chainsCommunity%>%
#   select(lp__,
#          #plot_mani intercepts (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.1.1, U.2.1.1, U.3.1.1, U.4.1.1,
#          U.1.2.1, U.2.2.1, U.3.2.1, U.4.2.1,
#          U.1.3.1, U.2.3.1, U.3.3.1, U.4.3.1,
#          U.1.4.1, U.2.4.1, U.3.4.1, U.4.4.1,
#          #plot_mani linear slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.1.2, U.2.1.2, U.3.1.2, U.4.1.2,
#          U.1.2.2, U.2.2.2, U.3.2.2, U.4.2.2,
#          U.1.3.2, U.2.3.2, U.3.3.2, U.4.3.2,
#          U.1.4.2, U.2.4.2, U.3.4.2, U.4.4.2,
#          #plot_mani quad slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#          U.1.1.3, U.2.1.3, U.3.1.3, U.4.1.3,
#          U.1.2.3, U.2.2.3, U.3.2.3, U.4.2.3,
#          U.1.3.3, U.2.3.3, U.3.3.3, U.4.3.3,
#          U.1.4.3, U.2.4.3, U.3.4.3, U.4.4.3,
#          #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
#          D.1.1.1, D.2.1.1, D.3.1.1, D.4.1.1,
#          D.1.1.2, D.2.1.2, D.3.1.2, D.4.1.2,
#          D.1.1.3, D.2.1.3, D.3.1.3, D.4.1.3,
#          #richness intercept, linear, and quad slopes (center digit): 2=richness
#          D.1.2.1, D.2.2.1, D.3.2.1, D.4.2.1,
#          D.1.2.2, D.2.2.2, D.3.2.2, D.4.2.2,
#          D.1.2.3, D.2.2.3, D.3.2.3, D.4.2.3,
#          #MAP intercept, linear, and quad slopes (center digit): 1=MAP
#          E.1.1.1, E.2.1.1, E.3.1.1, E.4.1.1,
#          E.1.1.2, E.2.1.2, E.3.1.2, E.4.1.2,
#          E.1.1.3, E.2.1.3, E.3.1.3, E.4.1.3,
#          #MAT intercept, linear, and quad slopes (center digit): 2=MAT
#          E.1.2.1, E.2.2.1, E.3.2.1, E.4.2.1,
#          E.1.2.2, E.2.2.2, E.3.2.2, E.4.2.2,
#          E.1.2.3, E.2.2.3, E.3.2.3, E.4.2.3,
#          #overall intercept, linear, and quad slopes
#          mu.1.1, mu.2.1, mu.3.1, mu.4.1,
#          mu.1.2, mu.2.2, mu.3.2, mu.4.2,
#          mu.1.3, mu.2.3, mu.3.3, mu.4.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:mu.4.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))
# 
# write.csv(chainsCommunity2, 'bayesian_output_summary_10 yr_03302017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_10 yr_03302017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,3640:5031]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 3640:5031])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,3640:5031]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 3640:5031])'] <- 'sd'
# 
# chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
#   #split names into parts
#   separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
#   select(-B)%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          id=as.integer(id))%>%
#   #if 95% confidence interval overlaps 0, then set mean to 0
#   mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
#   #spread by variable
#   select(variable, id, parameter, mean)%>%
#   spread(key=parameter, value=mean)
# 
# write.csv(chainsFinal, 'bayesian_output_mean sd_10 yr_03132017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_10 yr_03132017.csv')


# chainsByResource <- as.data.frame(chainsCommunity[,3640:5031])
# chainsByResource$rowID <- seq.int(nrow(chainsByResource))
# chainsByResource <- chainsByResource%>%
#   gather(key=parameter, value=value, 1:1392)%>%
#   separate(parameter, c('level', 'variable', 'study1', 'parameter'))%>%
#   mutate(study=as.integer(study1))%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')))%>%
#   select(-level)%>%
#   left_join(expInfo2, by='study')%>%
#   left_join(trtDetail, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
#   mutate(resource_mani=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
#                               ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
#                                      ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
#                                             ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
#                                                    ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'nuts:CO2',
#                                                           ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'nuts:dro',
#                                                                  ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'nuts:irr',
#                                                                         ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'CO2:dro',
#                                                                                ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'CO2:irr',
#                                                                                       ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'nuts:CO2:dro', 
#                                                                                              ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr', 'other'))))))))))))%>%
#   mutate(resource_mani_combo=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
#                                     ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
#                                            ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
#                                                   ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
#                                                          ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'multiple',
#                                                                 ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'multiple',
#                                                                        ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'multiple',
#                                                                               ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'multiple',
#                                                                                      ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'multiple',
#                                                                                             ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'multiple', 
#                                                                                                    ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'multiple','other'))))))))))))%>%
#   group_by(rowID, variable, parameter, resource_mani_combo)%>%
#   summarise(value_mean=mean(value))%>%
#   ungroup()%>%
#   group_by(variable, parameter, resource_mani_combo)%>%
#   summarise(value=mean(value_mean), sd=sd(value_mean))%>%
#   #if 95% confidence interval overlaps 0, then set mean to 0
#   mutate(lower=value-2*sd, upper=value+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, alt_mean=ifelse(diff==-2, 0, value))%>%
#   ungroup()
#   
# write.csv(chainsByResource, 'plot mani_equations_10 yr_by resource manipulated.csv')
chainsByResource <- read.csv('plot mani_equations_10 yr_by resource manipulated.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  select(-irrigation, -drought)%>%
  #get resource mani information
  left_join(trtDetail)%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*9+quadratic*9^2)*(0.2057881)+(0.3846015),
                    ifelse(variable=='dispersion', (intercept+linear*9+quadratic*9^2)*(0.102365)+(-0.008532193),
                           ifelse(variable=='evenness', (intercept+linear*9+quadratic*9^2)*(0.1133359)+(0.03247852), (intercept+linear*9+quadratic*9^2)*(0.2572354)+(-0.141815)))))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2057881)+(0.3846015),
                    ifelse(variable=='dispersion', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.102365)+(-0.008532193),
                           ifelse(variable=='evenness', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1133359)+(0.03247852), (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2572354)+(-0.141815)))))%>%
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
                                                                                                      ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'multiple','other'))))))))))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,', '*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=',
         curve6=') +',
         color=ifelse(resource_mani_combo=='nuts', '#00660044', ifelse(resource_mani_combo=='drought', '#FF660044', ifelse(resource_mani_combo=='irrigation', '#00009944', ifelse(resource_mani_combo=='CO2', '#66006644', ifelse(resource_mani_combo=='other', '#99999944', '#33333344'))))),
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_10 yr_by resource manipulated.csv', row.names=F)


###main figure
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,32,5), labels=seq(1,33,5)) +
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.4908421992764 + 0.20181995626785*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(-0.750291956005 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(-0.8518752203 + 0.297409347564*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(-1.01522138165 + 0.29310715573055*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.16509092053111*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0.6888899734575 + 0.20298956005215*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0.798510776725 + 0.2253241491793*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#99999944') +
  stat_function(fun=function(x){(0.5844242377036 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0.56152759687695 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0.715746645421 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(-1.0337314033 + 0.19202297259*x + -0.0114026560909*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(-0.89136532205 + 0.1515411133901*x + -0.0110741105499*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,15), colour='#99999944') +
  stat_function(fun=function(x){(-0.75879647041 + 0.188323053075*x + -0.01083990054825*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(-1.0273422061 + 0.114871597940115*x + -0.00828826965158*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(-1.1725787327 + 0.11851317969064*x + -0.00904815973439*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(-0.85214521615 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,13), colour='#33333344') +
  stat_function(fun=function(x){(-0.92621264645 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,13), colour='#00660044') +
  stat_function(fun=function(x){(-1.4439100762 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,13), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.11620932525*x + -0.0019268554205592*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0.4149312584789 + 0.168785579325*x + -0.0038144733955*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-0.47015705998675 + 0.33562228525*x + -0.0079827955985*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-1.10942363335 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-0.4107089902362 + 0.27307125575*x + -0.0065561864405*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-0.615923591265 + 0*x + 0.00193596974859045*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0.643063550535 + 0.07760922355745*x + -0.0028576254986825*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0.48771347633735 + 0.088764087239*x + -0.0019470628016205*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.0905507033283*x + -0.00463066030474*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0.3974811189326 + 0.1781327048*x + -0.0077005176595*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + -0.0028796913196725*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(-0.918611244545 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.439844810474655 + 0.177371681804604*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.66642233172 + 0.17552583931115*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.838861328365 + 0.1720236125108*x + -0.01738611015114*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.5015956388*x + -0.0265390328183*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.48457446660015 + 0.49531179745*x + -0.0274724715265*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-1.03827742945 + 0.144637196255935*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.6261907132645 + 0.33323563977*x + -0.0142304402551695*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.40950605335*x + -0.0168679113256*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.747156127965 + 0.1375623822615*x + -0.0099507595998455*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(-0.67263696637 + 0.128040647146*x + -0.008491321282235*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(-0.7435588945 + 0.1653522067775*x + -0.0097977127750085*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.83654252175 + 0.1560642793025*x + -0.009079220584*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.7056204456 + 0.127213995523326*x + -0.00892182686461*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.632720709075 + 0.13152360731955*x + -0.0086736986510035*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.62554288201 + 0.11639822436255*x + -0.00828336966821*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.7970959188 + 0.1640276008275*x + -0.009271148878465*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.8309470484 + 0.1806523688895*x + -0.010877400555425*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#99999944') +
  stat_function(fun=function(x){(-0.56323314140582 + 0.143330746506*x + -0.00932673336693*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(-0.84216516905 + 0.1112827340688*x + -0.00754763216884*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(-0.87676118275 + 0.10853069099508*x + -0.00777899655017*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(-0.82970204995 + 0.1590104887645*x + -0.0104240009223955*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(-0.72512985225 + 0.173464983025*x + -0.0108329814943*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(-0.85060427955 + 0.1041541556493*x + -0.0077086277604135*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,23), colour='#00660044') +
  stat_function(fun=function(x){(1.379505524 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,22), colour='#99999944') +
  stat_function(fun=function(x){(1.5628925102 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,22), colour='#00660044') +
  stat_function(fun=function(x){(-1.2360883103 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#99999944') +
  stat_function(fun=function(x){(-0.6078045360055 + 0.30995200829*x + -0.016313491232867*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.281494822656*x + -0.0147351820874575*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.1304968177917*x + -0.00512550413098125*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.10315899983684*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.119062533283766*x + -0.005285041289252*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0.13742229645675*x + -0.004524500069549*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.12424792667515*x + -0.004920566989038*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.10958260662308*x + -0.00627994569077*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.171976246675*x + -0.0079616737108*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0.09757902994569*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.45167528584195 + 0*x + -0.005228455784154*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.133244920290425*x + -0.00846010680925*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.11187847101335*x + -0.0059004921755705*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(-1.21719074195 + 0.123686410891105*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(-1.00657519415 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(-1.0711795326 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-1.14360617055 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-0.4595139019435 + 0.14770322052235*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-0.49797985395195 + 0.137266633183455*x + -0.0082129353077165*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-1.0503079193 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-0.8167313508 + 0.1730707615626*x + -0.008593963949661*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-0.6272967848035 + 0.2049414712905*x + -0.009753033677945*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(-1.2190415073 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(-1.01848607365 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(-1.0517646741 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(-0.69093262194 + -0.112896442702585*x + 0.0182891363734*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.41096822365728 + 0*x + 0.015088139836325*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.138752415393865*x + 0.019395718822*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.3864861205903 + 0*x + 0.017123835946815*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0.01660787539887*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.605456047247 + -0.13650302661256*x + 0.01861560303475*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.124234255142255*x + 0.01983735086*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.493537169763 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.262054806395*x + -0.0090756190144775*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.4652111984337 + 0.137049756800645*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.582182516235 + 0.120090699807133*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.16807270226575*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.34044412307365 + 0.21631478166985*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.12091624168163*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-1.2886400206 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-1.21644386225 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-1.2175558215 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-1.10809839105 + 0.40412509815*x + -0.0170539497810344*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-1.37724243225 + 0.2250126297865*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-1.4710675651 + 0.27528427403*x + -0.0132224269101945*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.5476199339684 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,7), colour='#66006644') +
  stat_function(fun=function(x){(-1.58801385725 + 0*x + -0.0151399369432*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-1.48562640075 + 0*x + -0.0156726950181256*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#66006644') +
  stat_function(fun=function(x){(-0.91117302275 + 0*x + -0.015500662049425*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,9), colour='#33333344') +
  stat_function(fun=function(x){(-1.0205881237 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(-0.45244530833715 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(-1.896516512 + 0.22871196834*x + -0.01241765277522*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(-0.600236475753 + 0*x + 0*x^2)*(0.2057881)+(0.3846015)}, size=0.5, xlim=c(0,8), colour='#00660044') +
  #last five are the main effect lines
  #nuts
  stat_function(fun=function(x){(-0.325356889 + 0.129974400*x + -0.004959197*x^2)*(0.2057881)+(0.3846015)}, size=3, xlim=c(0,31), colour='#006600') +
  #drought
  stat_function(fun=function(x){(-0.896263937 + 0.190173013*x + -0.011121278*x^2)*(0.2057881)+(0.3846015)}, size=3, xlim=c(0,16), colour='#FF6600') +
  #irrigation
  stat_function(fun=function(x){(-0.904014506 + 0.114073836*x + -0.006894904*x^2)*(0.2057881)+(0.3846015)}, size=3, xlim=c(0,19), colour='#000099') +
  #CO2
  stat_function(fun=function(x){(-1.033312431 + 0.079216706*x + -0.005920939*x^2)*(0.2057881)+(0.3846015)}, size=3, xlim=c(0,27), colour='#660066') +
  #other
  stat_function(fun=function(x){(-0.420719741 + 0.101626744*x + -0.006172811*x^2)*(0.2057881)+(0.3846015)}, size=3, xlim=c(0,25), colour='#999999') +
  #multiple
  stat_function(fun=function(x){(-0.783176507 + 0.122829034*x + -0.008837079*x^2)*(0.2057881)+(0.3846015)}, size=3, xlim=c(0,14), colour='#333333')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,32,5), labels=seq(1,33,5)) +
  # scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.4, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,15), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,13), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,13), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,13), colour='#66006644') +
  stat_function(fun=function(x){(-1.0227385438175 + 0.1048601278162*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-1.76988482405 + 0.125339227320795*x + -0.00355597664131445*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-2.142746722 + 0.1369520473526*x + -0.003825216924288*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0.99973667419255 + -0.1472159280145*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(-0.55785701193642 + 0.1786680046275*x + -0.0062924528258*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.270042260105*x + -0.014595911684*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.17836342140035*x + -0.00759385043143*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.13896153907506*x + -0.0059994521078965*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.124325962889565*x + -0.0058223291950165*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-1.32032316345 + -0.23384101857245*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.97386770801955 + -0.181521894470695*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.182153264337095*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.29570092276775*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.285477382252*x + -0.012926108423395*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,23), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.14136483258662*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,22), colour='#99999944') +
  stat_function(fun=function(x){(0 + -0.2564021315685*x + 0.008563905107391*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,22), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.80591145288065 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.159003391307295*x + 0.00818796983059015*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0.7713612198097 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.791480286515 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.8760410485065 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.639569397386935 + -0.161386510831565*x + 0.007517287843581*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0.7815273329672 + -0.152947404980009*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.213173930111341*x + 0.009472614348514*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.7844044143418 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.267397573485*x + 0.0124898367598*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0.7183089492562 + -0.19837153945075*x + 0.00855885622320145*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.7059874206941 + -0.18491114366655*x + 0.008591146762272*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.83916829085922 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.19681492758675*x + 0.008696389707866*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0.0094302477990865*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.174842630119305*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0.7755375466195 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.43543366995 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.49753834895 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.201182274455 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.45901171525 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.43784563035 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.5700634078 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.25095043590295*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.165486912469075*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2252198742292*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.1909247692295*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2755329522795*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.177226970358815*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,7), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#66006644') +
  stat_function(fun=function(x){(0.852875064391984 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,9), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0.20064788097831*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.102365)+(-0.008532193)}, size=0.5, xlim=c(0,8), colour='#00660044') +
  #last five are the main effect lines
  #nuts
  stat_function(fun=function(x){(0.230747108 + -0.044168548*x + 0.000000000*x^2)*(0.102365)+(-0.008532193)}, size=3, xlim=c(0,31), colour='#006600') +
  #drought
  stat_function(fun=function(x){(0.000000000 + 0.000000000*x + 0.000000000*x^2)*(0.102365)+(-0.008532193)}, size=3, xlim=c(0,16), colour='#FF6600') +
  #irrigation
  stat_function(fun=function(x){(0.08 + 0.000000000*x + 0.000000000*x^2)*(0.102365)+(-0.008532193)}, size=3, xlim=c(0,19), colour='#000099') +
  #CO2
  stat_function(fun=function(x){(-0.08 + 0.000000000*x + 0.000000000*x^2)*(0.102365)+(-0.008532193)}, size=3, xlim=c(0,27), colour='#660066') +
  #other
  stat_function(fun=function(x){(0.260889270 + -0.074217997*x + 0.000000000*x^2)*(0.102365)+(-0.008532193)}, size=3, xlim=c(0,25), colour='#999999') +
  #multiple
  stat_function(fun=function(x){(0.16 + 0.000000000*x + 0.000000000*x^2)*(0.102365)+(-0.008532193)}, size=3, xlim=c(0,14), colour='#333333')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,32,5), labels=seq(1,33,5)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.904070399025 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(1.141762047545 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0.9424828391365 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.19136259759875*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.226817187739935*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2317981506799*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#99999944') +
  stat_function(fun=function(x){(0 + -0.22814622864205*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2493060069023*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2553606594289*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0.7651286574865 + -0.13330278663625*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(0.7072293988405 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,15), colour='#99999944') +
  stat_function(fun=function(x){(0.58264169147545 + -0.205984464721*x + 0.01093571523446*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(0.86593406789 + -0.14475213239009*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(1.12055404945 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,13), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,13), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,13), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-0.88230663995 + -0.083548671907605*x + 0.0033646384206445*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-1.0982859689 + -0.086660795873*x + 0.003496843888997*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.25120350594*x + 0.00611707393983*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(1.17491924065 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2281158682*x + 0.00561191287*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0.5140719085576 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(-0.7730281589915 + -0.1154215094425*x + 0.003639976592675*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(-0.9199716254 + -0.1060765460057*x + 0.0028225482127635*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0.576441101008 + -0.098510984009*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.15287865127*x + 0.0042217926226055*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0.51503745975295 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0.6081186761615 + -0.22231659519005*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.844531742675 + -0.386918479405*x + 0.019448374828893*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.36758523383*x + 0.0184238093679093*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.80340252484465 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.3021899219335*x + 0.01551586450246*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.48035590054455 + -0.357087597285*x + 0.017601211026567*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.5152146894798 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2491667935175*x + 0.0156670827055845*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.8087868118995 + -0.287817421805*x + 0.017466514471745*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.770884125345 + 0*x + -0.0092037605676655*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(0.6396980654 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + -0.012332995140385*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0.137474418417945*x + -0.0126913217864265*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.60635747892685 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.6294717880405 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.793117230705 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.49393793246515 + 0*x + -0.01175958112026*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.42476432486825 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#99999944') +
  stat_function(fun=function(x){(0.5304507547179 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0.5416006695435 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0.142152180949275*x + -0.012340766156514*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.6716987350395 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0.50300826370405 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(0.5796753155075 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0.61043180001532 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,23), colour='#00660044') +
  stat_function(fun=function(x){(0.4064317942097 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,22), colour='#99999944') +
  stat_function(fun=function(x){(0 + -0.14222582022985*x + 0.005707978776321*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,22), colour='#00660044') +
  stat_function(fun=function(x){(0.6211145296205 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.55576984287379 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2333009336465*x + 0.01195302101295*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.738858110691 + -0.12024883844692*x + 0.00850681677963*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.1488763616535*x + 0.009165463554925*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0.466545478727938 + -0.17248122544196*x + 0.011509543374785*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.27399475557*x + 0.017921510932*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.18650461806445*x + 0.0121933956855*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + -0.1908930922385*x + 0.01161661468145*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.179797143699*x + 0.009465044274675*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0.520434596970435 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.14297911104549*x + 0.0060817665894535*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(0.515320448395 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0.625976870182 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0.511036714164989 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0.39510489742145 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0.48952214518965 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0.9992604928 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0.95237930815 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0.6848826995505 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.185911182266949*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.706548568622 + -0.1667383994838*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.59927944507995 + -0.142838041780565*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.283454502355*x + 0.0135987533174749*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.610623271344 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.232389942283*x + 0.012906361165232*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.943266665765 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.43935111972025 + -0.147189626499065*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.7937216622525 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.830373752459 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.5182198287428 + -0.167451831317035*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.842105918475 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.05823318275 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.6563789114184 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.602814567748845 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.7204629731815 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(1.09267912635 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.977933183945 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.820067057435 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0.864254861862 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,7), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0.944993249925 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#66006644') +
  stat_function(fun=function(x){(1.2756902505 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,9), colour='#33333344') +
  stat_function(fun=function(x){(0.6589411864895 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0.6376345197767 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0.47960377393075 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2572354)+(-0.141815)}, size=0.5, xlim=c(0,8), colour='#00660044') +
  #last five are the main effect lines
  #nuts
  stat_function(fun=function(x){(0.277233094 + -0.117613815*x + 0.003976229*x^2)*(0.2572354)+(-0.141815)}, size=3, xlim=c(0,31), colour='#006600') +
  #drought
  stat_function(fun=function(x){(0.673885174 + -0.169643626*x + 0.008320132*x^2)*(0.2572354)+(-0.141815)}, size=3, xlim=c(0,16), colour='#FF6600') +
  #irrigation
  stat_function(fun=function(x){(0.682671511 + 0.000000000*x + 0.000000000*x^2)*(0.2572354)+(-0.141815)}, size=3, xlim=c(0,19), colour='#000099') +
  #CO2
  stat_function(fun=function(x){(0.656242424 + 0.000000000*x + 0.000000000*x^2)*(0.2572354)+(-0.141815)}, size=3, xlim=c(0,27), colour='#660066') +
  #other
  stat_function(fun=function(x){(0.524630724 + -0.053112357*x + 0.000000000*x^2)*(0.2572354)+(-0.141815)}, size=3, xlim=c(0,25), colour='#999999') +
  #multiple
  stat_function(fun=function(x){(0.568822743 + 0.000000000*x + -0.006913742*x^2)*(0.2572354)+(-0.141815)}, size=3, xlim=c(0,14), colour='#333333')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,32,5), labels=seq(1,33,5)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00009944') +
  stat_function(fun=function(x){(-0.86834212636605 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,5), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.2088132728845*x + 0.01268623802172*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(0 + -0.1849716910602*x + 0.01243748134638*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,15), colour='#99999944') +
  stat_function(fun=function(x){(-0.915700454578 + -0.239322126755*x + 0.01372176212865*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,15), colour='#FF660044') +
  stat_function(fun=function(x){(0 + -0.188381965948*x + 0.012040131715885*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(0 + -0.157201021238102*x + 0.011166982087465*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,15), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,13), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,13), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,13), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0.5520020902329 + 0*x + -0.002953899220556*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0.9555735385785 + 0*x + -0.0034394107278565*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.17829799036*x + -0.003838238814585*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(-0.718784896975 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.192110103865*x + -0.0053096370722*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,30), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.10746614787425*x + -0.00385421599687*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.11866831496245*x + -0.004272835675444*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,29), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,26), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(1.7369346066 + 0.4857808613*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.370879356905*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.796438265852 + 0.1768399274822*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.21637574243895*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.32639305532*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.6030359416798 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.169249411248275*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.2462250984135*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,12), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0.004590938530526*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,23), colour='#00660044') +
  stat_function(fun=function(x){(0 + -0.10262474635277*x + 0.0047624894496081*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,22), colour='#99999944') +
  stat_function(fun=function(x){(0 + -0.10610281889391*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,22), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#99999944') +
  stat_function(fun=function(x){(-0.78183513215905 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.7859525279255 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,21), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,18), colour='#00009944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,11), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#99999944') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.7163774221098 + 0.152680815219192*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.5536894554103 + 0.185932797711425*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.8024939799287 + 0.152371795217828*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.59692767408795 + 0.176586132136795*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.858685571455 + 0.144051499359404*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.7037585754927 + 0.152008807833235*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.62018493885405 + 0.15805000502692*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.263177104562229*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.218506288120975*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(0 + 0.24699936080726*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,10), colour='#00660044') +
  stat_function(fun=function(x){(-0.73360504880456 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,7), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#00660044') +
  stat_function(fun=function(x){(-0.63896718045102 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,9), colour='#33333344') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0.7017074084306 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,16), colour='#66006644') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1133359)+(0.03247852)}, size=0.5, xlim=c(0,8), colour='#00660044') +
  #last five are the main effect lines
  #nuts
  stat_function(fun=function(x){(-0.161581820 + 0.073352163*x + 0.000000000*x^2)*(0.1133359)+(0.03247852)}, size=3, xlim=c(0,31), colour='#006600') +
  #drought
  stat_function(fun=function(x){(-0.678248060 + -0.224067700*x + 0.013204000*x^2)*(0.1133359)+(0.03247852)}, size=3, xlim=c(0,16), colour='#FF6600') +
  #irrigation
  stat_function(fun=function(x){(0.000000000 + 0.000000000*x + 0.000000000*x^2)*(0.1133359)+(0.03247852)}, size=3, xlim=c(0,19), colour='#000099') +
  #CO2
  stat_function(fun=function(x){(-0.267163164 + 0.000000000*x + 0.000000000*x^2)*(0.1133359)+(0.03247852)}, size=3, xlim=c(0,27), colour='#660066') +
  #other
  stat_function(fun=function(x){(-0.193784650 + 0.000000000*x + 0.000000000*x^2)*(0.1133359)+(0.03247852)}, size=3, xlim=c(0,25), colour='#999999') +
  #multiple
  stat_function(fun=function(x){(0.000000000 + 0.000000000*x + 0.000000000*x^2)*(0.1133359)+(0.03247852)}, size=3, xlim=c(0,14), colour='#333333')

# print(evennessPlot) #export at 1200x1000


#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000




###overall responses from bayesian output --------------------------------------------------------
###summary stats from bayesian output --------------------------------------------------------
# #gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(#plot_mani intercepts (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#     U.1.1.1, U.2.1.1, U.3.1.1, U.4.1.1,
#     U.1.2.1, U.2.2.1, U.3.2.1, U.4.2.1,
#     U.1.3.1, U.2.3.1, U.3.3.1, U.4.3.1,
#     U.1.4.1, U.2.4.1, U.3.4.1, U.4.4.1,
#     #plot_mani linear slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#     U.1.1.2, U.2.1.2, U.3.1.2, U.4.1.2,
#     U.1.2.2, U.2.2.2, U.3.2.2, U.4.2.2,
#     U.1.3.2, U.2.3.2, U.3.3.2, U.4.3.2,
#     U.1.4.2, U.2.4.2, U.3.4.2, U.4.4.2,
#     #plot_mani quad slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
#     U.1.1.3, U.2.1.3, U.3.1.3, U.4.1.3,
#     U.1.2.3, U.2.2.3, U.3.2.3, U.4.2.3,
#     U.1.3.3, U.2.3.3, U.3.3.3, U.4.3.3,
#     U.1.4.3, U.2.4.3, U.3.4.3, U.4.4.3,
#     #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
#     D.1.1.1, D.2.1.1, D.3.1.1, D.4.1.1,
#     D.1.1.2, D.2.1.2, D.3.1.2, D.4.1.2,
#     D.1.1.3, D.2.1.3, D.3.1.3, D.4.1.3,
#     #richness intercept, linear, and quad slopes (center digit): 2=richness
#     D.1.2.1, D.2.2.1, D.3.2.1, D.4.2.1,
#     D.1.2.2, D.2.2.2, D.3.2.2, D.4.2.2,
#     D.1.2.3, D.2.2.3, D.3.2.3, D.4.2.3,
#     #MAP intercept, linear, and quad slopes (center digit): 1=MAP
#     E.1.1.1, E.2.1.1, E.3.1.1, E.4.1.1,
#     E.1.1.2, E.2.1.2, E.3.1.2, E.4.1.2,
#     E.1.1.3, E.2.1.3, E.3.1.3, E.4.1.3,
#     #MAT intercept, linear, and quad slopes (center digit): 2=MAT
#     E.1.2.1, E.2.2.1, E.3.2.1, E.4.2.1,
#     E.1.2.2, E.2.2.2, E.3.2.2, E.4.2.2,
#     E.1.2.3, E.2.2.3, E.3.2.3, E.4.2.3,
#     #overall intercept, linear, and quad slopes
#     mu.1.1, mu.2.1, mu.3.1, mu.4.1,
#     mu.1.2, mu.2.2, mu.3.2, mu.4.2,
#     mu.1.3, mu.2.3, mu.3.3, mu.4.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:mu.4.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(CI=sd*2)%>%
#   separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
#   mutate(parameter=ifelse(level=='mu', predictor, parameter), predictor=ifelse(level=='mu', 'overall', predictor))%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          predictor=ifelse(level=='D'&predictor==1, 'ANPP', ifelse(level=='D'&predictor==2, 'rrich', ifelse(level=='E'&predictor==1, 'MAP', ifelse(level=='E'&predictor==2, 'MAT', ifelse(level=='U'&predictor==1, 'plot mani 2', ifelse(level=='U'&predictor==2, 'plot mani 3', ifelse(level=='U'&predictor==3, 'plot mani 4', ifelse(level=='U'&predictor==4, 'plot mani 5', 'overall')))))))))%>%
#   select(level, parameter, variable, predictor, predictor, median, sd, CI)
# 
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_10 yr_03302017.csv')

chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_10 yr_03302017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)


meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-1.1, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('\nMean Change') +
  annotate('text', x=3.45, y=-1.1, label='(a)', size=10, hjust='left')

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='dispersion' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.22, 0.4), breaks=seq(-0.2, 0.6, 0.2)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Dispersion\nChange') +
  annotate('text', x=3.45, y=-0.22, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='richness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.8), breaks=seq(-0.3, 0.5, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Richness\nChange') +
  annotate('text', x=3.45, y=-0.3, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.5, 0.25), breaks=seq(-0.3, 0.3, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.5, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 2400x500



#mean plots --------------------------------------------------------
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.3, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#dispersion plots --------------------------------------------------------
dispersionIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.3, 1.15) +
  coord_flip()

dispersionSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

dispersionQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#richness plots --------------------------------------------------------
richnessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.3, 1.15) +
  coord_flip()

richnessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

richnessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#evenness plots --------------------------------------------------------
evennessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.3, 1.15) +
  coord_flip()

evennessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

evennessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
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
#export at 2400x2000



###by resource mani - raw data--------------------------------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

rawTrt <- rawData%>%
  filter(experiment_length>9)%>%
  filter(treatment_year==experiment_length)%>%
  select(site_code, project_name, community_type, treatment, plot_mani, rrich, anpp, MAT, MAP, experiment_length, treatment_year, mean_change, dispersion_change, SimpEven_change, S_PC)%>%
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
                                                                                             ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr', 'other'))))))))))))%>%
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
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='mean_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.10), name='Mean Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(0, 0.6), xlim=c(1,6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.6, label='(a)', size=12, hjust='left') +
  annotate('text', x=1, y=0.54, label='ab*', size=10) +
  annotate('text', x=2, y=0.57, label='a*', size=10) +
  annotate('text', x=3, y=0.29, label='b*', size=10) +
  annotate('text', x=4, y=0.4, label='b*', size=10) +
  annotate('text', x=5, y=0.39, label='b*', size=10) +
  annotate('text', x=6, y=0.37, label='b*', size=10)

dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='dispersion_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Dispersion Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.1, 0.1), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.1, label='(b)', size=12, hjust='left') +
  annotate('text', x=1, y=-0.065, label='ab', size=10) +
  annotate('text', x=2, y=-0.056, label='ab', size=10) +
  annotate('text', x=3, y=-0.105, label='a*', size=10) +
  annotate('text', x=4, y=0.09, label='ab', size=10) +
  annotate('text', x=5, y=0.05, label='b', size=10) +
  annotate('text', x=6, y=-0.055, label='ab', size=10)

richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='S_PC', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.1), name='Richness Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.35, 0.17), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.17, label='(c)', size=12, hjust='left') +
  annotate('text', x=1, y=0.169, label='a', size=10) +
  annotate('text', x=2, y=-0.32, label='b*', size=10) +
  annotate('text', x=3, y=-0.16, label='a', size=10) +
  annotate('text', x=4, y=-0.29, label='ab*', size=10) +
  annotate('text', x=5, y=-0.2, label='b*', size=10) +
  annotate('text', x=6, y=-0.125, label='a', size=10)

evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='SimpEven_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Evenness Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.1, 0.15), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.15, label='(d)', size=12, hjust='left') +
  annotate('text', x=1, y=0.04, label='ab', size=10) +
  annotate('text', x=2, y=0.08, label='a*', size=10) +
  annotate('text', x=3, y=0.135, label='ab', size=10) +
  annotate('text', x=4, y=0.08, label='a', size=10) +
  annotate('text', x=5, y=-0.07, label='b*', size=10) +
  annotate('text', x=6, y=-0.07, label='b*', size=10)

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 2400 x 1600


###by magnitude of resource manipulated---------------------------------
#N addition
meanNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=mean_change)) +
  geom_point(size=5) +
  # scale_x_log10() +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=0.4, y=1, label='(a)', size=12, hjust='left')

dispersionNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=dispersion_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  # scale_x_log10() +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=0.4, y=0.4, label='(b)', size=12, hjust='left')

richnessNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=S_PC)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  # scale_x_log10() +
  scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=1.5, label='(c)', size=12, hjust='left')

evennessNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=SimpEven_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  # scale_x_log10() +
  scale_y_continuous(breaks=seq(-0.2, 0.6, 0.2), name='Evenness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.6, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600

#H2O change
meanPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=mean_change)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=-20, y=0.65, label='(a)', size=12, hjust='left')

dispersionPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=dispersion_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=-20, y=0.2, label='(b)', size=12, hjust='left')

richnessPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=S_PC)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks=seq(-1,2,0.25), name='Richness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-20, y=0.6, label='(c)', size=12, hjust='left')

evennessPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=SimpEven_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks=seq(-0.4, 0.5, 0.1), name='Evenness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-20, y=0.2, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600






#look at number replicates for dispersion results (all factors actually) -- doesn't make a difference --------------------------------------
reps <- read.csv('SpeciesRelativeAbundance_Dec2016.csv')%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year, plot_id)%>%
  summarise(mean=mean(relcov))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year)%>%
  summarise(rep_num=n())%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(rep_num=mean(rep_num))

dispersionReps <- rawTrt%>%
  left_join(reps, by=c('site_code', 'project_name', 'community_type', 'treatment'))

meanRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=mean_change)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Mean Change') +
  scale_x_continuous(breaks=seq(0,50,5)) +
  coord_cartesian(xlim=c(2,25)) +
  annotate('text', x=0, y=1, label='(a)', size=12, hjust='left')
dispersionRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=dispersion_change)) +
  geom_point() +
  geom_hline(yintercept=0) +
  xlab('Number of Relicates') +
  ylab('Dispersion Change') +
  scale_x_continuous(breaks=seq(0,50,5)) +
  coord_cartesian(xlim=c(2,25)) +
  annotate('text', x=0, y=0.35, label='(b)', size=12, hjust='left')
richnessRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=S_PC)) +
  geom_point() +
  geom_hline(yintercept=0) +
  xlab('Number of Relicates') +
  ylab('Richness Change') +
  scale_x_continuous(breaks=seq(0,50,5)) +
  coord_cartesian(xlim=c(2,25)) +
  annotate('text', x=0, y=1.4, label='(c)', size=12, hjust='left')
evennessRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=SimpEven_change)) +
  geom_point() +
  geom_hline(yintercept=0) +
  xlab('Number of Relicates') +
  ylab('Evenness Change') +
  scale_x_continuous(breaks=seq(0,50,5)) +
  coord_cartesian(xlim=c(2,25)) +
  annotate('text', x=0, y=0.6, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanRepPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionRepPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessRepPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessRepPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600







###look at five factor manipulations for mean change --------------------------------------------------------
#can't do this anymore because these experiments were not 10+ years (no NIN, no TRA, and no treatment with 4 factors including N)

