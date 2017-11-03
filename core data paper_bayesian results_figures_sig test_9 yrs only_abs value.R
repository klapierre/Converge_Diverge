library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

#kim's laptop
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34),
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
expRaw <- read.csv('ExperimentInformation_May2017.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('ForBayesianAnalysis_8yr_absvalue_May2017.csv')

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

#only run to generate initial chains files
#raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_abs_disp\\mv_abs_disp_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_abs_disp\\mv_abs_disp_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_abs_disp\\mv_abs_disp_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_abs_disp\\mv_abs_disp_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)


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
# write.csv(chainsCommunity2, 'bayesian_output_abs_summary_10142017.csv')

chainsCommunity2 <- read.csv('bayesian_output_abs_summary_10142017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments --------------------------------------------------------
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
# write.csv(chainsFinal, 'bayesian_output_abs_mean sd_10142017.csv')

chainsFinal <- read.csv('bayesian_output_abs_mean sd_10142017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,', '*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=grey) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_abs value.csv', row.names=F)



###main figure

# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
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
  stat_function(fun=function(x){(-0.4908744773237 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.662853612952915 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62560381011855 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.29954432648985*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.52687675906965 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.385463172135675*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.605219419543 + 0.3911639231555*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.044265054995 + 0.328622563760405*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.1605885158 + 0.40706304403762*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.60071579962115 + 0.3632807290597*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.2372340192 + 0.39969904044705*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.17468482075 + 0.4277216837*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.39560314905 + 0.40668569006245*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.61999086687665 + 0.28928178978205*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.866405986845 + 0.28042865077349*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7063094923995 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.61164676059835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.680513493822 + 0.347278623919265*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6134422268676 + 0.3911685152804*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.77837242516795 + 0.40550998380535*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.56280325435487 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3711562394016*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.34062873628975*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.303472351798105*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9141077327835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.7988085473995 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.1244230019 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.90187773412 + 0.2858146529815*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.799419520761 + 0.40317534392995*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.9078999366935 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6539032606858 + 0.3427953608136*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.87815574576525 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.85265623135 + 0.37340882076345*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.9414962331 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.81460894385 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.88418988465 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.86314064265 + 0.368034919077398*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.9758463119 + 0.30931736734345*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.58211354918742 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.55414678196915 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4953668063994 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5874612429283 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.42994229262445*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.567402891925 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.632764657291905 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.823275272685 + 0.40231923534555*x + -0.051203254467255*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8611418225255 + 0.4210492382325*x + -0.046138936319515*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.3568821759 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.27947017151435*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.8054204896*x + -0.074840247546*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.78977843555*x + -0.0792662761145*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.66448327940596 + 0.519927889574*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.2256904885 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.9891459282 + 0.69108851623*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.58553300465105 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.9671380102*x + -0.09343025643*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.9543557099*x + -0.086300068855*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6912812583*x + -0.036782635902637*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.84676550825*x + -0.049861316269925*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5340753785562 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.52649754699955 + 0.510762556348*x + -0.044730161057635*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.462777925575*x + -0.03557246203782*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.02962845405*x + -0.071216969907*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.73512203684295 + 1.1363109401*x + -0.0899966587885*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.610099224475*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.476221989932765 + 0.64385278388*x + -0.037802626899827*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.787459627518 + 0.317808057412245*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7075509181968 + 0.2993266923446*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.887346069405 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.840670515005 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.00057759754 + 0.33106684615985*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.8723997396 + 0.44842984339495*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.21183682 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.37704406795 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.5225681934 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.5145702875 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.45807869135 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.48538982785 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.4310404455 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.51297753225 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.35837726886545*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.9571421106028 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.6791844506497 + 0.4359230930985*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.1515242414 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.8023465252565 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.53164092447508 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.737258112385 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5807881117257 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.57763516834702 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7082592470845 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4933736010423 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.65416288982895 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.51053542097245 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.751152659715 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6215004880077 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.582797176819435 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62963710677295 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6348848333815 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.712635022621 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32437000521885*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.303528835282515*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5437469952755 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5440988209059 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.636696720841 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.60482919576295 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.90743103725 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-1.36791184725 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.18686629435 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.19457675765 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.209317875 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.12489623745 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.08752008525 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.2347393378 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1273800353 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.27161309065 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.585919752605665 + 0.28159922147925*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5772874963759 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.244085465998455*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5752778435097 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5061030455703 + 0.219319188695*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31022627697485*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5173165122867 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.51598164602645 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.4511418918163 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.217830268613905*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.632856734512 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32495264548066*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36463857455533*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.305399552329135*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.35448577005615*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.349624832928285*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5510975271226 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.321472245765*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.332825154950645*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(2.503710904 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(2.4488842995 + 0.31155640056855*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.656268969533 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.378179437412*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.476151642591*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.28133955682531*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3105719610398*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.9339752166655 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.801275937554 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8982370730685 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8154047585674 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8390178152005 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.866890511835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.84582335622 + 0.3439168897423*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.68189045615459 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7270760475165 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.77747397233085 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.13760843055 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0916364283 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.01540049995 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9176207090665 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.84507557849 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8500418681625 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.888539596405 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4800297228035*x + -0.0451800967082495*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4676843291195*x + -0.05155139087918*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8108638664293 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7272597256573 + 0.53662960546*x + -0.047149465530565*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.58473756329852 + 0.48353494839*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.1475587824 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.00419397735 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.07866028205 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.1266041142 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.113299586 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.91318534698 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.420706856925*x + -0.035586676700815*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2718488516341*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38377364564175*x + -0.032550813248036*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4334910971375*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32324454783775*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.53910007476*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.536778539992*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.86109750113175 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.348121097828575*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31210343932445*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.816322049303 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7744140734865 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.864429420415 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.692158938456835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7040590853309 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.4496329733649 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.54131718773625 + 0.27304219572415*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7407882700165 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6553707543895 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.55590947458495 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8737169608 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5757715963025 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7420964215 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.61955540574415 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7677226082075 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.633869242074 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.633799255097 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.49036983099485 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.13600131395 + 0.295809202213235*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.772180238552 + 0.4646014479085*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.963860917295 + 0.328646897906975*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.32691441422672*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.49498247919576 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.734113546723 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.322281738599235*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.72606436016035 + 0.45020147222285*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4277148510459*x + -0.06015954693225*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.67664636565315 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.63828057654655 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.63031003701755 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4792262529765*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.58567788695075 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3938701527981*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4451193236945*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4547546566325*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.800098121379 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.067867089563 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.91357457051 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.923843318145 + 0.472278653495*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.18150641645 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.15034150045 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.21706528745 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.1699568667 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.44634246045 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.8999197003725 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.9874125351255 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.924615078467 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.910901035015 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.927543174835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.9060883687665 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.91579729435 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.4497347178731 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.2950460362 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5503572634855 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.716288107939 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.724435237145 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.587911359598025 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.59099824713755 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5234696934282 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.44889905058652 + 0.31006046053175*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.740425354753795 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-1.4693986394 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.2904592609 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8389670963825 + 0.321924644231105*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.9296126501 + 0.3252826221609*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.901329994721 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.0469090949 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.045868279285 + 0.26126768507567*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.9212606197 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.75562274862835 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.91934007965 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.82821428514 + 0.3412385469865*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.97178305868 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.45771504543767*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.664047734941 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4086610395219*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4196066691096*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.543463068373975 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.68327695819405 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.6442707978528 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.268592960258365*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.33236912950305*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.58912330826044 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.546639374332*x + -0.06723130503848*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.73354897514945 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.65605104541105 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.271859721574874*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0403245368096295*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.322264492326164*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6194709138412 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31198865459052*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5231376369563 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1417218)+(0.2924074)}, size=0.5, xlim=c(0,2), colour='grey') +
  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  stat_function(fun=function(x){(-0.46141200 + 0.19014700*x + -0.01285265*x^2)*(0.1417218)+(0.2924074)}, size=3, xlim=c(0,8), colour='black')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(0,0.45))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.45, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.68197633233785 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.67560635464431 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.593811195761215 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6535292046165 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.65423634461545 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.57612994825965 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.641412403271955 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.67507074602 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7312541675985 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.68441625195515 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.62330566187828 + 0.40207240421075*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.327392810425 + 0.5748631388705*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.054790250893555*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.42234457924375*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.48747175363995*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38632537650475*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.1097384186975 + 0.47381074920715*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.34591944944385*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.37384772180916*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.88228500007925 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.86239119951195 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.75278584695 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7228395191904 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8958897347178 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.23781076858 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9695114055395 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.0729253983745 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.181250281115 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.34911346665 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.86866829108725 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.390014374595495*x + -0.064257034282*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4689607403887*x + -0.0719365599939*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.43892375525225*x + -0.07844919049125*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6544802432979 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.595588103184135 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6551281146279 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5812313422173 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.597969111004045 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4706396857164*x + -0.08714219072421*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.73511793346505 + 0.7695433088*x + -0.10176377759*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.698072210813249 + 0.5468731515244*x + -0.109290601345*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.51195628158735*x + -0.09641394833*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.6702240156475 + 0.633986585232*x + -0.105775265659907*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.6505529259145 + 0.663588012885*x + -0.10181754507*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7819890865555 + 0.58943657495*x + -0.109317035175*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.60487871660665 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5737877381749 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.74069457023925 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.685876623982677 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.754969773242115 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.048051232410555*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.05785120518959*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.71960270680595 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.49301070605045*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38268949093914*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.400868477253415*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.364148574698712*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.996830980143 + 0.492831131619*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.0918523396645 + 0.55372414501675*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.43971864471185*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3739484034871*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.41542447482415*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36692326284679*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6963484029293 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.60970926806055 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.354053734332035*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.361614731955619*x + -0.057263574408595*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3389995374615*x + -0.06066862157077*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5453662800385*x + -0.065993767576505*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.317126326541*x + -0.05943641691109*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3733983471865*x + -0.052072786849116*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.74471827057585 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.67562567250435 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.644003519913485 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.63154191120688 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.60318959479375 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.660241517195587 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.8085625720269 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0573189048880135*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.35537594901565*x + -0.057451085325901*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.060410042121485*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.70090720277775 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.580225237559545 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=0.5, xlim=c(0,2), colour='grey') +
  #estimated as mean across treatment lines
 stat_function(fun=function(x){(-0.24528650 + (0.0603264+0.33940800)*x + 0*x^2)*(0.05918496)+(0.06230041)}, size=3, xlim=c(0,8), colour='black')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(0,0.4))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change') +
  annotate('text', x=0, y=0.4, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.540661596961805 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.523122240974235 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5333944204834 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.501221325886825 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.527053374974812*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.63147477718295 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.511163396678755 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6491832198845 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.56587955103919 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62396243986795 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.557517717841 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.72630754464455 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.60368175759519 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.638101017394425 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6252059237379 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5514060715251 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.55243878620548 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.567596983350521 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.54587527364427 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.3880696472*x + -0.12318409742*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.7379170564*x + -0.1603499381*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.44758136751125*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5784855399419*x + -0.07799202725875*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6198229234659*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.24157875175*x + -0.115208958995*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.35528291345*x + -0.135069975025*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.62731652914185 + 0.70977730848665*x + -0.09726772282815*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.15011350495*x + -0.136539144015*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.39158121685*x + -0.122280765365*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.8618767835*x + -0.1858660221*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.661464996220915 + 1.25000814265*x + -0.1102843734935*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.4953566118*x + -0.1321096877775*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.96823907085*x + -0.08923703314813*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 1.16210457025*x + -0.117948639261*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7303415356465 + 0.51436601580435*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.52383677372133 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.56124236956165 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.63279370469575 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.663640938141885 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.5809337124985 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.6943549311774 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.5962440844923 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.57911470227885 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5216109316598 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.55386917994765 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.52752892775877 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.5165240702247 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.55180097069945 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.65641610291825 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.55535512522925 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5141078312468 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6112927381807 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5036381806099 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.54934120322335 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.577153406294481 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.57608774384035 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5426450386422 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.47585292804755 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.4549983115298 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.475124161029615 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.477813249969*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6017071838882 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.759841312349 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.58594244240578 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.510158024881634 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.51986766941197 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6761318570366 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.58875863493765 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5092174785157 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.40113455661385*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.631780880435865 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.64213614286735 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.55437451966415*x + -0.077802532759355*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.689180752657*x + -0.099794545124*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4859715504645*x + -0.0553541616193661*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.44589313535645*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5694206672236 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.57503069005975 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9370422954975 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.0315058810142 + 0.80262045375485*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.130215856604 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.6745599500089 + 1.11310021285*x + -0.1019054845138*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.43675007012936*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.47475385855*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.5258740356047 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.5054441042418 + 0.4232612519102*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36977478263303*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.60912243033965 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.57260620965155 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.57191630368565 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.641843287890446 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.5147405565217 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.50408459784875 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.58554088407195 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6721240715355*x + -0.068133170650674*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5150049333769*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.362959107747931*x + -0.0566116632801975*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.545761782105*x + -0.0606952294872*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.591057894399958 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.5413991031823 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 2.0411776675*x + -0.23414022116*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 1.01734549645*x + -0.129250033559895*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7497087430718 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.601418676081675 + 0*x + 0*x^2)*(0.1539214)+(0.1548579)}, size=0.5, xlim=c(0,2), colour='grey') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  stat_function(fun=function(x){(-0.23596100 + 0*x + 0*x^2)*0.1539214 + 0.1548579}, size=3, xlim=c(0,8), colour='black')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(0,1.3))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=1.3, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5798762407597 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.58817503818155 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.75592722824*x + -0.06689058990961*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.44992072817095 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.59317761708115 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.38185003034168 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.40961987159345 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.369358457323565 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4036825542284 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4751754413642 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.48168005954491 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.47401972813578 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.48597828235085 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.3942045082829 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5691852066867 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.57681645649665 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5507383552181 + 0.52666936518184*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6890138814666 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8669818344595 + 1.1808560144*x + -0.09751110262255*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.851341519027 + 2.0502855455*x + -0.15733421842*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5115492153547*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5623947556295 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.539839035111877 + 0.668348078401*x + -0.06908739650487*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.471510874306204*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7242048141597 + 1.32712692665*x + -0.103362272166*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 3.0036363365*x + -0.2903758213*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.532273047985815 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.70200103266885 + 0.556212669635*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.718449056936 + 0.95840031315*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.565642654665299*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5903151524733 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5993025650071 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.58167147636076 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7666794165495 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.625008962328 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.770129264145 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.6315712215049 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.77229289237 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.61754969242 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.768464596305 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.34905175648755 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.36019049508425 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.3709716924365 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.357195263070295 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.36268573994625 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.3526500183703 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.3718686740977 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.338416877271 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.504925959691365 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.47544703436605 + 0.46354336048995*x + -0.066184187283685*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.491615077568505 + 0.6356685700665*x + -0.0562628523384715*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.49770245563922 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-0.483099298199565 + 0.385382519985525*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.34933964748475 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.426830135478535 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.45975289810985 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.45236839156 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.444008141604355 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.367611208020015 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.4599418420718 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.367149939957195 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.34650233337595 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5614295907559*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.464932888345647 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.468217425507035 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.42339722230805*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.44210282974125 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.4505170349989 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.44518409497902 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.45091924930693 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.46866565043185 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.470685691697183 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.45013352271925 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.36420629159335*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.46973702790441 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.56683252055005 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.691248730715505 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.6031210444 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.55960464105 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.62724993405 + -0.59901347145915*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.55263924983355*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 1.03146947659*x + -0.09349339729504*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 1.1123870673*x + -0.095661849765175*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 1.52177706705*x + -0.1306830151744*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.660900545562*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6594701811826 + 0.627040913421*x + -0.06523759364021*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.45072295560985 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.45188180123535 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.60232207129155 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.4535241932762 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6053114573418 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.595622368864 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.51339429093311 + 0.52243549671335*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.52844742192895 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.531631428760345 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.069598863802085*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.4560685946469 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5010532329614 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.504358693329 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.49117100111146 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5886948284657 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + -0.46920794535415*x + 0.06987892272951*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.528730821547111 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.4541187465629 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6267667737383*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.46039643309165 + 0.746400991738*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.45397772002265 + 0.60490234963965*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.47698264317955 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.47717457676175 + 0.6993050563494*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4228293351483 + 0.911697519997*x + -0.091307107674445*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.45470309383735 + 0.570635566461715*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4864353993845 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.37995660954155 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.39221506772338 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.4137132150038 + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=0.5, xlim=c(0,2), colour='grey') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  stat_function(fun=function(x){((-0.20869150) + 0*x + 0*x^2)*(0.06881917)+(0.06727672)}, size=3, xlim=c(0,8), colour='black')

# print(evennessPlot) #export at 1200x1000


#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000



# ##summary stats from bayesian output --------------------------------------------------------
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
# write.csv(chainsCommunitySummary, 'bayesian_output_abs value_summary_final plots_10142017.csv')

chainsCommunitySummary <- read.csv('bayesian_output_abs value_summary_final plots_10142017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)


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
  ylim(-1.2, 1.4) +
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
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.16,0.1)) +
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
  ylim(-1.2, 1.4) +
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
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.16,0.1)) +
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
  ylim(-1.2, 1.4) +
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
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.16,0.1)) +
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
  ylim(-1.2, 1.4) +
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
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.16,0.1)) +
  coord_flip()

#plot all together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(4,3))) 
print(meanIntPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(meanSlopePlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(meanQuadPlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 3))
print(dispersionIntPlot, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
print(dispersionSlopePlot, vp=viewport(layout.pos.row = 4, layout.pos.col = 2))
print(dispersionQuadPlot, vp=viewport(layout.pos.row = 4, layout.pos.col = 3))
print(richnessIntPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(richnessSlopePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessQuadPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(evennessIntPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessSlopePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(evennessQuadPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
#export at 2500x2000


#overall responses --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.8, 0.5, 0.4)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('\nMean Change') +
  annotate('text', x=3.45, y=-0.8, label='(a)', size=10, hjust='left')

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='dispersion' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.8, 0.5, 0.4)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Dispersion\nChange') +
  annotate('text', x=3.45, y=-0.8, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='richness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.8, 0.5, 0.4)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Proportion\nRichness Change') +
  annotate('text', x=3.45, y=-0.8, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.5), breaks=seq(-0.8, 0.5, 0.4)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.8, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 2400x500