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
expRaw <- read.csv('ExperimentInformation_May2017.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('ForBayesianAnalysis_8yr_May2017.csv')

rawData2<- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_disp=mean(dispersion_change), std_disp=sd(dispersion_change), mean_rich=mean(S_PC), std_rich=sd(S_PC), mean_even=mean(SimpEven_change), std_even=sd(SimpEven_change)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani))

#for table of experiment summarizing various factors
rawDataAll <- read.csv('ForBayesianAnalysis_May2017.csv')
expInfoSummary <- rawDataAll%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
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
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_3.csv', comment.char='#')
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
# write.csv(chainsCommunity2, 'bayesian_output_summary_03132017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_03132017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,12172:17235]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 12172:17235])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,12172:17235]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 12172:17235])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_03132017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_03132017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=8, 7, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*7+quadratic*7^2)*(0.1417218)+(0.2924074),
                    ifelse(variable=='dispersion', (intercept+linear*7+quadratic*7^2)*(0.0859427)+(0.0001271936),
                           ifelse(variable=='evenness', (intercept+linear*7+quadratic*7^2)*(0.09476219)+(0.01686738), (intercept+linear*7+quadratic*7^2)*(0.2127569)+(-0.04917495)))))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1417218)+(0.2924074),
                    ifelse(variable=='dispersion', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.0859427)+(0.0001271936),
                           ifelse(variable=='evenness', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.09476219)+(0.01686738), (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2127569)+(-0.04917495)))))%>%
  mutate(color=ifelse(variable=='mean'&20>=rrich, '#000AE544', ifelse(variable=='mean'&rrich>=21&40>=rrich, '#1F09C144', ifelse(variable=='mean'&rrich>=41&60>=rrich, '#3F079D44', ifelse(variable=='mean'&rrich>=61&80>=rrich, '#5F057944', ifelse(variable=='mean'&rrich>=81&100>=rrich, '#7F035544', ifelse(variable=='mean'&rrich>=101&120>=rrich, '#9F013144', ifelse(variable=='mean'&rrich>120, '#BF000D44', ifelse(variable=='dispersion'&-10>MAT, '#000BE544', ifelse(variable=='dispersion'&MAT>=-10&-6>=MAT, '#1B09C644', ifelse(variable=='dispersion'&MAT>=-5&-1>=MAT, '#3607A744', ifelse(variable=='dispersion'&MAT>=0&4>=MAT, '#51068844', ifelse(variable=='dispersion'&MAT>=5&9>=MAT, '#6D046944', ifelse(variable=='dispersion'&MAT>=10&14>=MAT, '#88034A44', ifelse(variable=='dispersion'&MAT>=15&19>=MAT, '#A3012B44', ifelse(variable=='dispersion'&MAT>=20, '#BF000D44', ifelse(variable=='richness'&-10>MAT, '#000BE544', ifelse(variable=='richness'&MAT>=-10&-6>=MAT, '#1B09C644', ifelse(variable=='richness'&MAT>=-5&-1>=MAT, '#3607A744', ifelse(variable=='richness'&MAT>=0&4>=MAT, '#51068844', ifelse(variable=='richness'&MAT>=5&9>=MAT, '#6D046944', ifelse(variable=='richness'&MAT>=10&14>=MAT, '#88034A44', ifelse(variable=='richness'&MAT>=15&19>=MAT, '#A3012B44', ifelse(variable=='richness'&MAT>=20, '#BF000D44', ifelse(variable=='evenness'&MAP<200, '#000BE544', ifelse(variable=='evenness'&MAP>=200&399>=MAP, '#1B09C644', ifelse(variable=='evenness'&MAP>=400&599>=MAP, '#3607A744', ifelse(variable=='evenness'&MAP>=600&799>=MAP, '#51068844', ifelse(variable=='evenness'&MAP>=800&999>=MAP, '#6D046944', ifelse(variable=='evenness'&MAP>=1000&1199>=MAP, '#88034A44', ifelse(variable=='evenness'&MAP>=1200&1399>=MAP, '#A3012B44', '#BF000D44')))))))))))))))))))))))))))))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,', '*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=', curve6=') +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_env varis.csv', row.names=F)




###main figure (Figure 1)
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1), labels=seq(1,8,1)) +
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.64024733638295 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.488681930397755 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.63562716630475 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.29782563845203*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.508896469869775 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.415248847614345*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6041614440291 + 0.3822881809366*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.0404712696 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.15198119535 + 0.41825092448825*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.5930329175765 + 0.4056781996953*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.2063201668 + 0.415507855447305*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.159185262 + 0.45300475445555*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.36440440055 + 0.42294590267225*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.62855618349158 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.8396302450125 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.6543847476725 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.6128743872171 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6898939695685 + 0.35429197593575*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.616335628058 + 0.4236649162232*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7895120887855 + 0.4423858937603*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.55836882539355 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38488511869354*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.90258034767 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.7969040973435 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.12470489455 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.883437257335 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.793052847063 + 0.4225796308355*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.905719516663 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.70446055167 + 0.4107247050225*x + -0.04277275760029*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.84144007566 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.8580373251 + 0.374682463318145*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.9487679388 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.8716149452 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.9601112235 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.8554466863 + 0.4092278625772*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-2.007535691 + 0.3251461551109*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.562828260667185 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5215806701598 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.56067794172115 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.511868075316755 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.426785648941789*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4308836995891*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.54683315228635 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.65843257450115 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.874473132289 + 0.464460788056*x + -0.0606729726642*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.919130657105 + 0.4893222453765*x + -0.05609443275928*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.33997440205 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.834089805*x + -0.07904263403*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.863717437*x + -0.0884692014335*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.595436920577 + 0.5404526583525*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.28448467515 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8512402394933 + 0.545252136497*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.546360039858685 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 1.0469351421*x + -0.10522885087*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.98268796055*x + -0.090281317024*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.60292642615*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.8482246418*x + -0.05148593338286*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.57652818033501 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.66072708359515 + 0.591811730825*x + -0.054114452278415*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5075916861105*x + -0.044538640351894*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.9764374002*x + -0.063353197519245*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.660618654434 + 1.18720517025*x + -0.0930828027895*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.68939629085*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5215966766385 + 0.6302661431*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.78362138923015 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.681340920811655 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.87490429491785 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.80783631181 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.05741269235 + 0.34286259542999*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.93806238582 + 0.486491337153*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.1822153664 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.4194157243 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.5374272127 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.52449751505 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.46990473395 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.48651868015 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.409678243 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.53825063205 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3415255205006*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.9847057342345 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7236987365215 + 0.44048342948535*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.2054938864 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.815325621835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.476901660083875 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7939772198314 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5534807398165 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.57714358929217 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6965931498598 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6700102934297 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.47519271091155 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.709949840866 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5951223870172 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.63960827353395 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5955250383177 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.630084853676 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.73765801202 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.39006615159525*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5950702149557 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.6732887365006 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.665908445337 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.90464180787 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.35347628875 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.18666244385 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1703816542 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.27774997245 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1438157459 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0905970146 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.19930473355 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1008417328 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.2327439058 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.56301503774365 + 0.27635257772542*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5858300087575 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.53337573131205 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.51214123711475 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31648580041321*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5630273162776 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5518576683515 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.486224077781895 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.640829963263 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.519688446330815 + 0.37395312682141*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.50635529153304 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.353696717333545*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(2.4585126085 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(2.4289245765 + 0.3450853578041*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.67622567923565 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38684734168305*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5258746851765*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.34052281735765*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.914550473344 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7098930725065 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8758408599235 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.55114598103535 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7735408924955 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8396771645685 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.82353957356615 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.787091566239 + 0.39695551318441*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.657989328348 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.784345000364 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7889710437075 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1541626133 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.10069410205 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.02396089032 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.924422140965 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8523810556855 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.899267315835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.914538909435 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5775587900992*x + -0.05798555101985*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.53753846594*x + -0.0601460470869*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.79888918318 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7771356977485 + 0.558870206935*x + -0.04978001069265*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6939797023205 + 0.5852201597575*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.1828296779 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.95437283815 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.04438116455 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.0881625908 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.13853570205 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.92415691181 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.47020508477805*x + -0.042529083257405*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3990032937115*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.41997797218585*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3439790870056*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.56873303273*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5196691133637*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.863084139137 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.345053431010085*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.77075338667515 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7262848774543 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.820755396414 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7401089263 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6160005894165 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7434826718355 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6191767349305 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.54866910975865 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.9099354402 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.554865714943545 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7749631804429 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.66270874012845 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6861354211825 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.673920553018 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5870653821672 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.10669866606 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.780822178831 + 0.480433168157*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7759721375725 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.7301447735725 + 0.587898241865*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.388277333520635*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.63040893739465 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6709646488215 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.64146926063565 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.539883573735*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6209767280741 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.43709852218035*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.488231498655*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.52827209535905*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.834056717924 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.145235207775 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.9853329519915 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.9689616196 + 0.51915260344705*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.254852126 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.2101929374 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.255184882 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.2176960643 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.48754552875 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.8928378935385 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.96908764767 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.96128336051895 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.94281305043 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.995801679865 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.97074233862 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.90017237231 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.31108000405 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.598554218421565 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7626513285855 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7323595441719 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.63878365238465 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.605357464215335 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.513509599426191 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.4781808396584 + 0.33448510021815*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.782482166055 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.3442609318 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.29163776645 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8994028785305 + 0.4918839384297*x + -0.0686632657554*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.99552738115 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.890253518005 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.01223698467 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.0197513733 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.953901131475 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7220492912805 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.9681154548 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.817704152895 + 0.3674053782481*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7875600266987 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.462309660937847*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6284755223422 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.397689731708415*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4139807023404*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.600706374890875 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.6988535437305 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.6584267857074 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.351590148223*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.534408018058*x + -0.065494113107935*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7757629449545 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.616602029279076 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.645848686319 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5337635087633 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.46463900 + 0.16878700*x + -0.01217110*x^2)*(0.1417218)+(0.2924074)}, size=3, xlim=c(0,7), colour='black')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1), labels=seq(1,8,1)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.4, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6621738204827 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.5554277454867 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.77513572578655 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6405538843545 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.68085765510305 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.896289388949 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.759074154044515 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.97521824701555 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.090764565217 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.71756591865 + -0.52048774086945*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.570885379086*x + -0.068200419436695*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.061852240175795*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0649001266357*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8186014597281 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.3770118666 + -0.46228141384804*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3909495558536*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5869933204394*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7375385885041 + -0.43400661439685*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4832758726024*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.563499000808*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.83536680034561 + -0.40471549052476*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9547670936875 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9454732802855 + -0.3863838263321*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.0730237991115 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.8311431190685 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.77218987835425 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7002949080861 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.63878460314168 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.436340160396965*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.607533366318992 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8780353312685 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.70399447406785 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.6346184762 + 0.8576772304032*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.2865780555 + 0.78810978752439*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.7342950917 + 0.769187711359*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.61471848965 + 0.8705471615685*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.32481412525 + 0.70589780070375*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.40072513205 + 0.983845260568*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.96536294935 + 0.7889409297335*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.278837525315 + 0.953546621522*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.77872497565 + 0.87219844500415*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.301289277 + 0.73487677364285*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.112357603989 + 1.0438964399941*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.78341953666385 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.644610651018298 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.649034529839 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7025078696719 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.876085854105305 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8676427394303 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8861745537036 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.84152747449893 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.81325991348575 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.82405593115712 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7354480661343 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.68337379748785 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.20396073223 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.059227775573405*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.09835908925 + 0.714481323335*x + -0.0971021287635*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.055621386165 + 0.5256977470424*x + -0.0986437440875*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7927984025184 + 0.466596668123296*x + -0.08258887134625*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.997100367906 + 0.5522317737544*x + -0.0870609867904*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.06000133979 + 0.6408527250085*x + -0.0994482877345*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.12537549463 + 0.5502877004489*x + -0.096818108675*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.052351304767625*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38794221899575*x + -0.0812352126738*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0621706359125115*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05972044207376*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05829219462349*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.72455263143885 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.87333449923505 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.39523688445095*x + 0.08126920225644*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0875491008913*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.56777543271395 + -0.6079495018675*x + 0.0835474836856*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.0122999054575 + -0.5481326519905*x + 0.058879548917385*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.721069183745*x + 0.0843268468735*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.155528957435 + -0.42641900831455*x + 0.06738874084706*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.601239727567805 + -0.6200862499415*x + 0.072203157383505*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.78437986699743 + -0.49229774377575*x + 0.0590873665208625*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6685691327386 + -0.59112893055*x + 0.090949707478*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.82440917176575 + -0.4418480936009*x + 0.0658822546074*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6038681918735*x + 0.07806538422958*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6932913141535 + -0.581233865*x + 0.07517537013848*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.773338508728785 + -0.55937768954*x + 0.07329419262645*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.76216321645091 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.8557822140398 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.89482266052015 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.23866389789 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.77460308427265 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.28050908265 + 0.4586510187349*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.5517737223 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.837435596631 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.672168785894045 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.955823686719 + 0.55772152702935*x + -0.073245233737055*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.370020264037905*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.753503869750145 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.86406337548339 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.70522294111075 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.91700081080135 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.90479129602375 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.840699253474 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1520368719325 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.53660413566365*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0184784502995 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4927452010747*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.9774768711205 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6827536514738 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7593988367617 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8868922112825 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.0924310293 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.67562189499579 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.287622851235 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.843903176639 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.929179726298 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.720071021049375 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.822476993245015 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.62576922995 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.4078570585 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.900991993762729 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.173152436175 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62764115344945 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4452723052879*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7532462266132 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4961572738725*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  #estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)}, size=3, xlim=c(0,7), colour='black')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1), labels=seq(1,8,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.589389680743665 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.62232032556459 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.68245355239295 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.782460125191 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.69266721835475 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.940316077675 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0689581032034375*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.07411319814631*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7329819647077 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.65600711679305 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.61466843972735 + 0*x + 0.0660937494489095*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.74977177763045 + 0*x + 0.062931826214152*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.57385127121627 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.614433375375 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.69266736531655 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.754551080148 + -0.4334889360784*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5552906740584*x + 0.061844958840813*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.737680288063685 + -0.39504162433626*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6687789431936 + -0.46399165680885*x + 0.057003812371945*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.1533093604665 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.70754057138943 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.966386661375 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.0366359640645 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7497038368324 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7710317842322 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.102345846685 + -0.812775976725*x + 0.08054312745365*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.79901527026985 + -1.46692366485*x + 0.134801578245*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.55280736725*x + 0.15128175203*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.99894963785*x + 0.09576012712*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.5668033389 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.6517501365 + -1.1263797554*x + 0.093231879764*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.12111268245*x + 0.1127409126195*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.15742004485*x + 0.1189935416125*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4439662211145*x + 0.0779380334605285*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9418862738215 + -0.94685514465*x + 0.11616703709*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.2217696644*x + 0.14662147177*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5369326410838*x + 0.056057137165815*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.28866151*x + 0.118082486685*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.4119806793*x + 0.13848017437*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.60047691675525 + -0.429373627009895*x + 0.0542192852502*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.12798460545*x + 0.1009776849095*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.38823523755*x + 0.12869390204*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6699175057516 + -0.9665484207*x + 0.108362728722*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -1.0279833269*x + 0.10802671099*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.57897377023855*x + 0.070902983084225*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.61093802311709 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.7170503182424 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(1.21658109445 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.662494400823635 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.59181610569265 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.19251152415 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6622798969089 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.42548992395 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.17072976898 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.975852326116 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.55752067115 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.738223402904 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.37420596825 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.955639863065 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.50535242455 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.901267004765 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.733597061141 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.57824563808945*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6173150670388*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4156272274258*x + -0.0759668125797*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.40607633839185*x + -0.05744560201276*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0527942389797295*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.052674698787755*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.0969934133115 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.933930958075 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.52819718305*x + -0.07875686787265*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.485980680519*x + -0.0674017962044*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32560952138114*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3702889431692*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3938323703972*x + -0.058699530576885*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4545524435242*x + -0.0686902426584*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.39310657095805*x + -0.050215642986315*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.339541941141282*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4560220034345*x + -0.058164245841615*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.313824119706155*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.913670678320445*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.669538767895972*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.8949262143955*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.711767537444109*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.90923710949205*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.69328364514485*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.80382244433225*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.71460743959585*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.48031453452587*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3845141571415*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.62830296192985 + -0.48154625574495*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 1.56588862457*x + -0.44187784991575*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 1.4821771444135*x + -0.427791369682*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.86489016710048 + 1.0417117131705*x + -0.3818221763965*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.69191714519394*x + -0.08792727304641*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9457243573795 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.715755004717*x + -0.091922310713315*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.78386575398195 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6945273838665 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.81284325951945 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8372457939915 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.050672533686905*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.387524858042225*x + -0.0898687368089*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0494132367127167*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05261120505018*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.387128912928505*x + -0.0810079406283*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.411899885702315*x + -0.0829399670597*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.099979736965 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.259849815125 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.17210582299 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.26990042195 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.012758200829 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.93084459008905 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.585223284047*x + 0.08987043848855*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.67439850107855 + -0.63110243652919*x + 0.09214411464975*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0732706709263*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.339052928534275*x + 0.0677028910302*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.37184495137305*x + 0.07020488486524*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6483735363735*x + 0.105700745115*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0487605560912696*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.51098844956175*x + 0.089176389327*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6672741995*x + 0.117373951865*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4606981759111*x + 0.097497096255*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.356095055706075*x + 0.059585593676635*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4015571893391*x + 0.0858844984415*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0775681537947*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.64403647755271 + -0.65874104603105*x + 0.0982600427191*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.936040988455 + -0.535284205823345*x + 0.0868996435306*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.507230354041*x + 0.081306949948675*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.71163972270915 + 0*x + 0.07731934411848*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.684932608285*x + 0.11013707268631*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.694255164820245*x + 0.11567138096235*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.09160644210101*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.69444967901667 + 0*x + 0.0920514567902575*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.900160317008 + 0*x + 0.0840314931337073*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.49127445116825*x + 0.086989117246795*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.6683912896575 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.65998982228715 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.4987918201 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.8327863498 + -2.9596040967*x + 1.05059510385*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(2.4350865435 + -2.1490210837165*x + 0.96653145695*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(2.43825012535 + -2.23552796207*x + 0.97500034885*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.52747329154695*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.580548854086*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.89235866654*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.57159639805645*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.54654344817915*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.8226631846285 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.1324287762 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.46295979891549*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.450546213086465*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.42149593727153*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.9015500451775 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.8868851251975 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.571762150243555 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.68409310318801 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.743049026976245 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.57305295151862 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.585182163171365 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.6426518814166 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.774237451788 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.7709412217016 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.93127921590405 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.63722872859605 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5985246107375*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.8633436278628 + 0.7018950657235*x + -0.1023564420365*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5736673027145*x + -0.08608559400675*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7209464214722 + 0.617215999408*x + -0.0834112565705*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9606071137405 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7399886490415 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.967382614935 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4681327689117*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.568060739471445 + -0.409347939124705*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.79704855596985 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.5905475876159 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.73115207599405 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.7744339799353 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.180182935147025*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.7795209679739 + -0.5703729335262*x + 0.18052359593959*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6910777804472*x + 0.19972138385785*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.116572577931207*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4912157762284*x + -0.1360216201439*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.75867203114615*x + -0.161716596478*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.68809647731 + 0.5415342468522*x + -0.143774566000105*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.8921617620686 + 0*x + -0.1358419567221*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.12582549598502*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.12551455423455*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.1224165681562*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.961487278 + 0.954884885035*x + -0.2040048165765*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.81597904682185*x + -0.1673161021772*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.126583661024959*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.51203557270055*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.73484587632255 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.51358662265 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9931105801945 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.33740820575 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.29502450 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,7), colour='black')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,7), breaks=seq(0,7,1), labels=seq(1,8,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.64906671069285*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.66030847451615 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6556026352075 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.56485449237735 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6594451394724 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.528748761550645 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5385836788896 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.61682836093035*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.93350124633365 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.953196056613 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.19482343327 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.902119997414 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.0521716083135 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.400861620731779*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.531672972858*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.30279256639 + 0.91501554379*x + -0.1069587027791*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.4513895664 + 1.46216797125*x + -0.14137296808*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.220852717145 + 1.959642172*x + -0.160603309485*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.44886152302205*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7202260872854 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.58212245559185*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6198337965671*x + -0.095886091188*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6495411340478 + 1.3904794865*x + -0.133732974853*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 2.1193586035*x + -0.20296584706*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.75605869008206 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.9995659399105 + 0.75761805139*x + -0.0622024013392535*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8602811280013 + 1.01700626767*x + -0.063089898324497*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5593587054111*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.7814292601614 + 0.47738936388704*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.74079538020005 + 0.500622542404005*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.46777631404658 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.47329184992015 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.46795680994235 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.559944836841 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4845160572172 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.42308755059015*x + -0.0781116955984*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4639299210367*x + -0.083726764635253*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.520595845847*x + -0.0898144715786*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.70480147541*x + -0.0845829027952*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.566058029713*x + -0.0861485259205*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.519391036762*x + -0.0928410129815*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.7787986077295 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62193643511075 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8767732011935 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6956057764083 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8345141816675 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7180845184296 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.86273816015 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.74139595430445 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.52541431265445 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5358349783497 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.669362791986 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5918492309136 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6392044787986 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.77537001517098 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.72776177494079 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.619771669669365*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.47444015687955*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.57339046203475 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.632113243580105 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5116087577352 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.58454384064525 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.54805535204115 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.59925362058945 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.637456772545055 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.652267785796*x + -0.09277436198435*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.52466172436537*x + -0.09396923974981*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6919373380965*x + -0.09929452809875*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.476239852572955 + 0*x + -0.0524981111393145*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.474608149227056 + 0.365177167287945*x + -0.06875109136559*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.058572373596785*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5098274702395*x + -0.080415748792*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.55857554994815 + 0.3610767696821*x + -0.061534454951435*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.534628832713365*x + -0.0806373364132*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.37194594380655*x + -0.06274479154785*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.558341892320055 + 0.4323289046062*x + -0.0734595934743655*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3721462653705*x + -0.05957674208923*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.37089455401725*x + -0.06232688823925*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.56082547629*x + -0.0841741930285*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7185264079282 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.029765876989 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.92602509709689 + -0.773963997471*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.75977711325065*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.69105570573505*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6909902289927*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.8245713895875*x + 0.083829166254727*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.9496242148358*x + 0.092921151893225*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.712180158669815 + -0.9126669783915*x + 0.092660253843125*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -1.16943191436*x + 0.125110893556165*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -1.2120039088*x + 0.12751219489727*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.71537188626655 + -1.2094400399*x + 0.12997866751755*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.49650887142987*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-2.1016904474 + 0.644322039805*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-2.0863835019 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-2.15313895785 + 0.751962763824*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.56660289276379 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.6298852649594 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.56058883024645 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.93029295492*x + -0.08848251822093*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.93304685226*x + -0.082124011054365*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 1.31300288995*x + -0.11003611868269*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.52422133025783*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.57732641850455*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.641339813777815*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.73992341814313*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6817494638708*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6161131162948*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6010977581895*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.39409080031798*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3755516577795*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4109655706004*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.54292430850915*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.60275199457195 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.56904865194735 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.769031518366015 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.71755262540125 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.67584463058336 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.754031931671 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.90098410469678 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8963583325625 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.67942814160285 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5239696937945*x + -0.0663240974723*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.40644109411975*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5696569277365 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.41714119428085*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.59070749323265 + 0.4948267821194*x + -0.063995186227825*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.65988159050645 + 0*x + -0.06470826575195*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62461722487393 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.664894475684975 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.444719745206362*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.37221068102603*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.6047470950353 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.7484065047074*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6145263025442 + 0.6135001051598*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.49381641077076*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.51046611514015*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='grey') +
  #overall mean line
  stat_function(fun=function(x){(-0.19162650 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,7), colour='black') +

# print(evennessPlot) #export at 1200x1000


#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000


###summary stats from bayesian output --------------------------------------------------------
# # gather summary stats needed and relabel them
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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_09122017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_09122017.csv')

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
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.15, 1.15) +
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
  ylim(-1.15, 1.15) +
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
  ylim(-1.15, 1.15) +
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
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.14,0.1)) +
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
  ylim(-1.15, 1.15) +
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
  scale_y_continuous(name='Mean Change') +
  stat_function(fun=function(x){(0.2003873 + 0.003866423*x)}, size=5) +
  xlab('') +
  annotate('text', x=0.4, y=0.7, label='(a)', size=12, hjust='left')

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
  scale_y_continuous(name='Dispersion Change') +
  xlab('') +
  annotate('text', x=0.4, y=0.4, label='(b)', size=12, hjust='left')


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
  scale_y_continuous(name='Richness Change') +
  stat_function(fun=function(x){(0.04821532 - 0.0001258736*1000 - 0.01420553*x + 0.00001198415*1000*x)}, size=5, color='#4793CF') +
  stat_function(fun=function(x){(0.04821532 - 0.0001258736*600 - 0.01420553*x + 0.00001198415*600*x)}, size=5, color='#2D5E88') +
  stat_function(fun=function(x){(0.04821532 - 0.0001258736*200 - 0.01420553*x + 0.00001198415*200*x)}, size=5, color='#153049') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.7, label='(c)', size=12, hjust='left') +
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
  scale_y_continuous(name='Evenness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.6, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
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
  scale_y_continuous(name='Mean Change') +
  stat_function(fun=function(x){(0.1820251 + 0.0002544999*x)}, size=5) +
  xlab('') +
  annotate('text', x=-80, y=0.6, label='(a)', size=12, hjust='left')

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
  xlab('') +
  annotate('text', x=-80, y=0.3, label='(b)', size=12, hjust='left')

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
  scale_y_continuous(name='Richness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.6, label='(c)', size=12, hjust='left')

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
  scale_y_continuous(name='Evenness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.3, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanH2OPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionH2OPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessH2OPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessH2OPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600






###comparing different resource manipulation types (Figure 3)
#get means and sd to backtransform
trtInteractionsData <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\treatment interactions_09072017.csv')%>%
  summarize(mean_mean_change=mean(mean_change), sd_mean_change=sd(mean_change), mean_dispersion_change=mean(dispersion_change), sd_dispersion_change=sd(dispersion_change), mean_richness_change=mean(S_PC), sd_richness_change=sd(S_PC), mean_evenness_change=mean(SimpEven_change), sd_evenness_change=sd(SimpEven_change))


#mean change
trtmean <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_mean_posteriors.csv', comment.char='#')*0.1471785+0.3233751
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
  scale_y_continuous(name='Mean Change') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr'), labels=c(expression(paste(CO[2], '(5)')), 'Drought (13)', 'Irrigation (24)', 'N (58)', 'P (16)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (27)', 'N*P*K (18)', 'N*P*K*Irr (2)')) +
  annotate('text', x=0.5, y=0.6, label='(a)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#dispersion
trtdispersion <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_disp_posteriors.csv', comment.char='#')*0.09374439+0.003188805
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
  scale_y_continuous(name='Dispersion Change') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr'), labels=c(expression(paste(CO[2], '(5)')), 'Drought (13)', 'Irrigation (24)', 'N (58)', 'P (16)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (27)', 'N*P*K (18)', 'N*P*K*Irr (2)')) +
  annotate('text', x=0.5, y=0.2, label='(b)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#richness
trtrichness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_richness_posteriors.csv', comment.char='#')*0.2397557-0.04606815
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
  scale_y_continuous(name='Richness Change') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr'), labels=c(expression(paste(CO[2], '(5)')), 'Drought (13)', 'Irrigation (24)', 'N (58)', 'P (16)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (27)', 'N*P*K (18)', 'N*P*K*Irr (2)')) +
  annotate('text', x=0.5, y=0.5, label='(c)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#evenness
trtevenness <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\manipulation\\trt_simp_posteriors.csv', comment.char='#')*0.1003621+0.01021106
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
  scale_y_continuous(name='Evenness Change') +
  xlab('') +
  scale_x_discrete(limits=c('CO2', 'Drought', 'Irrigation', 'N', 'P', 'CO2...N', 'N...Drought', 'N...Irr', 'N...P', 'N...P...K', 'N...P...K...Irr'), labels=c(expression(paste(CO[2], '(5)')), 'Drought (13)', 'Irrigation (24)', 'N (58)', 'P (16)', expression(paste(CO[2],'*N (3)')), 'N*Dro (2)', 'N*Irr (8)', 'N*P (27)', 'N*P*K (18)', 'N*P*K*Irr (2)')) +
  annotate('text', x=0.5, y=0.15, label='(d)', size=12, hjust='left') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

pushViewport(viewport(layout=grid.layout(2,2)))
print(meantrtPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersiontrtPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennesstrtPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnesstrtPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
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
  scale_y_continuous(breaks=seq(0, 0.7, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,0.7)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.37, label='a*', size=10) +
  annotate('text', x=2, y=0.52, label='b*', size=10) +
  annotate('text', x=3, y=0.73, label='c*', size=10)
#export at 800x800