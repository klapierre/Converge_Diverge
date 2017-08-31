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
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_3.csv', comment.char='#')
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
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*9+quadratic*9^2)*(0.1482058)+(0.298937),
                    ifelse(variable=='dispersion', (intercept+linear*9+quadratic*9^2)*(0.08790049)+(-0.000407613),
                           ifelse(variable=='evenness', (intercept+linear*9+quadratic*9^2)*(0.09597775)+(0.01730774), (intercept+linear*9+quadratic*9^2)*(0.2154907)+(-0.0546117)))))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1482058)+(0.298937),
                    ifelse(variable=='dispersion', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.08790049)+(-0.000407613),
                           ifelse(variable=='evenness', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.09597775)+(0.01730774), (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2154907)+(-0.0546117)))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,', '*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations.csv', row.names=F)



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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4880356285543 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.45923792482907 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.66586951009855 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5250315682174 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.640773827848 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.25757595905441*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5404936284174 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.354105382177295*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.58501463563885 + 0.329385026931662*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.946833895925 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.06753258455 + 0.39171663612874*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.53274747115145 + 0.37066027177706*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.11293496237 + 0.390476509535275*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.07359375015 + 0.41897335039325*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.2690057616 + 0.3842415707243*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.57174617140715 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.7622962892045 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.556564261254215 + 0.284172173943291*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5130751539688 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.61133178370005 + 0.34179419479215*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.54961286543425 + 0.39357808503216*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.715750952466 + 0.4085068781549*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.48884344518602 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.350421070533528*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.879202746665 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.788128931895 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1126880988 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.857621485935 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.745880234975 + 0.3723046539623*x + -0.035862083629651*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7115185410245 + 0.3426156965076*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.92519866867 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.88570267216 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.7986788608 + 0.365683846496*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.89930168615 + 0.29888422502044*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.82709309685 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.9124644872 + 0.295353690383739*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.80059329745 + 0.38452493035505*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.9521091336 + 0.32936794020945*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5664881354532 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.54468584106293 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.57590565859974 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5224090584855 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.33094150653223*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.63724653995945 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3446931623096*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6772028994697 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.879267719765 + 0.426533186987*x + -0.0567289584015*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.891482077605 + 0.4376403480415*x + -0.05019061284904*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.32672833505 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.6900315236*x + -0.05430989986385*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.7085753277*x + -0.060783108118*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.62676064123035 + 0.56005171707*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.40115740505 + 0.321279605582*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.8921372899 + 0.55004182298*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.60442181683515 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.88780323865*x + -0.078102983715*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.83805710255*x + -0.0656276586365*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.507059363135*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.8736890528*x + -0.0657927888493*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.40027870857725*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.649230614237 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.46565402556075*x + -0.035498743594585*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.69331081598125 + 0.53222225511*x + -0.043198227624305*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.86655768075*x + -0.04909573545405*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.6720643846135 + 1.03928743295*x + -0.073473626097*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.72097740255675 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.5809344828565 + 0.6439503155*x + -0.04045133113151*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.82754417995*x + -0.05036201631015*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.66841962832905 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8527217965805 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8073073483255 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.02749081015 + 0.3114764691818*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.907599806565 + 0.4381684525315*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.18579151305 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-1.4873675445 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.50798589985 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.4985000855 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.44363512285 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.44607135655 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.38744870575 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.3663485167 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.33914335261305*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.991657156804964 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-1.18297494855 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.82647590537 + 0.47308191158356*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.61342820406575 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.71828055368295 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6523232776005 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6340397307809 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.64646618814935 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7678890957375 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.44928618360565 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6915362501855 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5037623508437 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.73332944547425 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.624235437207 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.819492130605 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5741683339895 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.835275162505 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.50426867556035 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.33137245443635*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.614248607641 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.516730983633655 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7392555661608 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7460876879173 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.888737356736 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.33598065275 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.17165028345 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.15810118085 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.2619205229 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.12286318245 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.07017430955 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.19095790475 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0816227472 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.22274499225 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.408328092359 + 0.231755048385745*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.470930527484245 + 0.24196338316523*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.612870220387365 + 0.3651126941225*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6484944537313 + 0.28246777617299*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.561282114364435 + 0.288705174617875*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.46660340846955 + 0.3215538056104*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.404829261695226 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5971892894991 + 0.2425671891591*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.54773065242509 + 0.28148699924265*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.42379447802098 + 0.38589992999935*x + -0.032878481926839*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.617633919675 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.63396078365845 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5298283600336 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.41023536606585 + 0.275587890491925*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.68807784945 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.462996537590114 + 0.33801590382765*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.32069122238825*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.33452968771336*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.4423361684107 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.32704649944815*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(2.3043407545 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(2.2814483665 + 0.3038814420699*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7631329567425 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3697228966433*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.49763273079*x + -0.0323267902299335*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.292370794742792*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.92306848161 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.72393072243445 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.884606580725 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.5736004394843 + 0.42944034441591*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4391409288354*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.38880574130675*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.46195347854826*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.38626585371208*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.752442794255 + 0.3647692460401*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.660562541379125 + 0.44316535631445*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.83678489877 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8371405065135 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7715725830195 + 0.35573242984766*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.658072700062 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.79126127892265 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7839301764775 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1622831834 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1091236255 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0178561054 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.935032467515 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8549774165095 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.907660894145 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.93421619196 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.488630364945*x + -0.044090580393505*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.46224586233*x + -0.04755410511955*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.82382181399 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7603169171735 + 0.5125529340375*x + -0.04540333621505*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6492559946635 + 0.498493412892*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.2129809612 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.97448486735 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.052444138935 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.06792693082 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1545712664 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.92449460768 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3124968912317*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.30719715125555*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.378067457439*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.28027838951668*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.51201993208305*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.479289454639*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8981385871295 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.320496420727425*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7765072122655 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.734830033571 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.831037327493 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7550802613665 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.64102682612735 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.48708753321275 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.41973714407126 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.772574465535 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.65105313611095 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.418890359659795 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.60204901108655 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.95212735247 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.602831105542 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.807787795125 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.605733988339 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.601441928654 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.5841201303643 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.492602842541835 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.35946934170765*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.32334554936632*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.358740164401275*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8588641425264 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.52750478771749 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5616235748196 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4882346838266 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.347539340497245*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.45807813085835*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.61909555000965 + 0.528252358188*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.57806552093895 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.47496311067171 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6284117768875 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6506795971163 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6430504667765 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4578904943705*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-0.614159056602 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.405884988828*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.43269221980695*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4597967411885*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.46708207478092 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(-0.46728964538461 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-0.853492667715 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.12260238275 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.953420847978 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.966985848825 + 0.52334233593*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.2945552981 + 0.29210399815135*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.25230898529 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.3187370982 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.265550022 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.4366574081 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.887154655295 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-0.97236633026 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.951578012435 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9515059461825 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-1.011186813335 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.022892044955 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8837789814 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.29293042665 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5925101211796 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.751202547325 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.727521622149 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.64731408511235 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.608845764287 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.551844762709765 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.478897367635985 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.47290897499684 + 0.29215034033605*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.781779264362 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.35181957385 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.28341122845 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.90582354604 + 0.430196197072*x + -0.0587985084865*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9663360588 + 0.2776868759236*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8910040071285 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00408007915 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.983147238475 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.96284991331 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7422240275285 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.985229257345 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.80612636933 + 0.3188153821085*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8277391998985 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.425065480698705*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6664371596737 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3763263829166*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3847670102878*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6240408041048 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6319372429875 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.59199667485231 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.323991115433855*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.510842046086*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.829413326517 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6187163333037 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.45005513657574 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.719255739718 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5512219744646 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1482058)+(0.298937)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines (plot mani 1-4 staggered by intercept so lines don't overlap)
  #mani1
  stat_function(fun=function(x){(-0.47095450 + 0.16878700*x + -0.01154170*x^2)*(0.1482058)+(0.298937)}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.49 + 0.16878700*x + -0.01154170*x^2)*(0.1482058)+(0.298937)}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.51 + 0.16878700*x + -0.01154170*x^2)*(0.1482058)+(0.298937)}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.53 + 0.16878700*x + -0.01154170*x^2)*(0.1482058)+(0.298937)}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.47095450 + (0.16878700+0.45703950)*x + (-0.01154170-0.03508225)*x^2)*(0.1482058)+(0.298937)}, size=3, xlim=c(0,8), colour='#EC1804')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.97521824701555 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.090764565217 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.71756591865 + -0.52048774086945*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.570885379086*x + -0.068200419436695*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.061852240175795*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0649001266357*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8186014597281 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.3770118666 + -0.46228141384804*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3909495558536*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5869933204394*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(-0.78341953666385 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.644610651018298 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.649034529839 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7025078696719 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.876085854105305 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8676427394303 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8861745537036 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.84152747449893 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.81325991348575 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.82405593115712 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7354480661343 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.68337379748785 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.20396073223 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.059227775573405*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.09835908925 + 0.714481323335*x + -0.0971021287635*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.055621386165 + 0.5256977470424*x + -0.0986437440875*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7927984025184 + 0.466596668123296*x + -0.08258887134625*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.997100367906 + 0.5522317737544*x + -0.0870609867904*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.06000133979 + 0.6408527250085*x + -0.0994482877345*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.12537549463 + 0.5502877004489*x + -0.096818108675*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.052351304767625*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38794221899575*x + -0.0812352126738*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0621706359125115*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05972044207376*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05829219462349*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + -0.39523688445095*x + 0.08126920225644*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0875491008913*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.56777543271395 + -0.6079495018675*x + 0.0835474836856*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.0122999054575 + -0.5481326519905*x + 0.058879548917385*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.721069183745*x + 0.0843268468735*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.155528957435 + -0.42641900831455*x + 0.06738874084706*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.601239727567805 + -0.6200862499415*x + 0.072203157383505*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.78437986699743 + -0.49229774377575*x + 0.0590873665208625*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.6685691327386 + -0.59112893055*x + 0.090949707478*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.82440917176575 + -0.4418480936009*x + 0.0658822546074*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6038681918735*x + 0.07806538422958*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.6932913141535 + -0.581233865*x + 0.07517537013848*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.773338508728785 + -0.55937768954*x + 0.07329419262645*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.840699253474 + 0*x + 0*x^2)*(0.0859427)+(0.0001271936)}, size=0.5, xlim=c(0,8), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)}, size=3, xlim=c(0,8), colour='black')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.589389680743665 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.62232032556459 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.68245355239295 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.782460125191 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.69266721835475 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.940316077675 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0689581032034375*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.07411319814631*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7329819647077 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.65600711679305 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.61466843972735 + 0*x + 0.0660937494489095*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.74977177763045 + 0*x + 0.062931826214152*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.57385127121627 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.614433375375 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.69266736531655 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.754551080148 + -0.4334889360784*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5552906740584*x + 0.061844958840813*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.737680288063685 + -0.39504162433626*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6687789431936 + -0.46399165680885*x + 0.057003812371945*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.1533093604665 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.70754057138943 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.966386661375 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.0366359640645 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.7497038368324 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7710317842322 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(1.102345846685 + -0.812775976725*x + 0.08054312745365*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.79901527026985 + -1.46692366485*x + 0.134801578245*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.55280736725*x + 0.15128175203*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.99894963785*x + 0.09576012712*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.5668033389 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(1.6517501365 + -1.1263797554*x + 0.093231879764*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.12111268245*x + 0.1127409126195*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.15742004485*x + 0.1189935416125*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.4439662211145*x + 0.0779380334605285*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.9418862738215 + -0.94685514465*x + 0.11616703709*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.2217696644*x + 0.14662147177*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.5369326410838*x + 0.056057137165815*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.28866151*x + 0.118082486685*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.4119806793*x + 0.13848017437*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.60047691675525 + -0.429373627009895*x + 0.0542192852502*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.12798460545*x + 0.1009776849095*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.38823523755*x + 0.12869390204*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.6699175057516 + -0.9665484207*x + 0.108362728722*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.0279833269*x + 0.10802671099*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.57897377023855*x + 0.070902983084225*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.61093802311709 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.7170503182424 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.21658109445 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.662494400823635 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0.59181610569265 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.19251152415 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6622798969089 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.42548992395 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.17072976898 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.975852326116 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.55752067115 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.738223402904 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.37420596825 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.955639863065 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.50535242455 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.901267004765 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.733597061141 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.57824563808945*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.6173150670388*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4156272274258*x + -0.0759668125797*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.40607633839185*x + -0.05744560201276*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.0527942389797295*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.052674698787755*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(1.0969934133115 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.933930958075 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.52819718305*x + -0.07875686787265*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.485980680519*x + -0.0674017962044*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.32560952138114*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3702889431692*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3938323703972*x + -0.058699530576885*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4545524435242*x + -0.0686902426584*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.39310657095805*x + -0.050215642986315*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.339541941141282*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4560220034345*x + -0.058164245841615*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.313824119706155*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.913670678320445*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.669538767895972*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.8949262143955*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.711767537444109*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.90923710949205*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.69328364514485*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.80382244433225*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.71460743959585*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.48031453452587*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3845141571415*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0.62830296192985 + -0.48154625574495*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.56588862457*x + -0.44187784991575*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 1.4821771444135*x + -0.427791369682*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.86489016710048 + 1.0417117131705*x + -0.3818221763965*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.69191714519394*x + -0.08792727304641*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9457243573795 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.715755004717*x + -0.091922310713315*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.78386575398195 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6945273838665 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.81284325951945 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.8372457939915 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.050672533686905*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.387524858042225*x + -0.0898687368089*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0494132367127167*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.05261120505018*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.387128912928505*x + -0.0810079406283*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.411899885702315*x + -0.0829399670597*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(1.099979736965 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.259849815125 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.17210582299 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.26990042195 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.012758200829 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.93084459008905 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.585223284047*x + 0.08987043848855*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.67439850107855 + -0.63110243652919*x + 0.09214411464975*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0732706709263*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.339052928534275*x + 0.0677028910302*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.37184495137305*x + 0.07020488486524*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.6483735363735*x + 0.105700745115*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0487605560912696*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.51098844956175*x + 0.089176389327*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.6672741995*x + 0.117373951865*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.4606981759111*x + 0.097497096255*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.356095055706075*x + 0.059585593676635*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4015571893391*x + 0.0858844984415*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0775681537947*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.64403647755271 + -0.65874104603105*x + 0.0982600427191*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.936040988455 + -0.535284205823345*x + 0.0868996435306*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.507230354041*x + 0.081306949948675*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.71163972270915 + 0*x + 0.07731934411848*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.684932608285*x + 0.11013707268631*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.694255164820245*x + 0.11567138096235*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.09160644210101*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.69444967901667 + 0*x + 0.0920514567902575*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.900160317008 + 0*x + 0.0840314931337073*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.49127445116825*x + 0.086989117246795*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6683912896575 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.65998982228715 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.4987918201 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.8327863498 + -2.9596040967*x + 1.05059510385*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(2.4350865435 + -2.1490210837165*x + 0.96653145695*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(2.43825012535 + -2.23552796207*x + 0.97500034885*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.52747329154695*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.580548854086*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.89235866654*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.57159639805645*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.54654344817915*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.8226631846285 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.1324287762 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.46295979891549*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.450546213086465*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.42149593727153*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.9015500451775 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.8868851251975 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.571762150243555 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0.68409310318801 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.743049026976245 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0.57305295151862 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.585182163171365 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6426518814166 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.774237451788 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7709412217016 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.93127921590405 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.63722872859605 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5985246107375*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.8633436278628 + 0.7018950657235*x + -0.1023564420365*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5736673027145*x + -0.08608559400675*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.7209464214722 + 0.617215999408*x + -0.0834112565705*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.9606071137405 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7399886490415 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.967382614935 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4681327689117*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.568060739471445 + -0.409347939124705*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.79704855596985 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5905475876159 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.73115207599405 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.7744339799353 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.180182935147025*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.7795209679739 + -0.5703729335262*x + 0.18052359593959*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.6910777804472*x + 0.19972138385785*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.116572577931207*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4912157762284*x + -0.1360216201439*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.75867203114615*x + -0.161716596478*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.68809647731 + 0.5415342468522*x + -0.143774566000105*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8921617620686 + 0*x + -0.1358419567221*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.12582549598502*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.12551455423455*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.1224165681562*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.961487278 + 0.954884885035*x + -0.2040048165765*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.81597904682185*x + -0.1673161021772*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.126583661024959*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.51203557270055*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.73484587632255 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.51358662265 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9931105801945 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.33740820575 + 0*x + 0*x^2)*(0.2127569)+(-0.04917495)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #plot mani 1-4 staggered to prevent overlap
  #mani1
  stat_function(fun=function(x){(0.29502450 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(0.32 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(0.34 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(0.36 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(0.38 + (-0.0622582-0.52537800)*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#EC1804')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.64906671069285*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.66030847451615 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6556026352075 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.56485449237735 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6594451394724 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.528748761550645 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.5385836788896 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.61682836093035*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.93350124633365 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.953196056613 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.19482343327 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.902119997414 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0521716083135 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.400861620731779*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.531672972858*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.30279256639 + 0.91501554379*x + -0.1069587027791*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.4513895664 + 1.46216797125*x + -0.14137296808*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-1.220852717145 + 1.959642172*x + -0.160603309485*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.44886152302205*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7202260872854 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.58212245559185*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.6198337965671*x + -0.095886091188*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.6495411340478 + 1.3904794865*x + -0.133732974853*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 2.1193586035*x + -0.20296584706*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.75605869008206 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.9995659399105 + 0.75761805139*x + -0.0622024013392535*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.8602811280013 + 1.01700626767*x + -0.063089898324497*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.5593587054111*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7814292601614 + 0.47738936388704*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.74079538020005 + 0.500622542404005*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.46777631404658 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.47329184992015 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.46795680994235 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.559944836841 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4845160572172 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.42308755059015*x + -0.0781116955984*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4639299210367*x + -0.083726764635253*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.520595845847*x + -0.0898144715786*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.70480147541*x + -0.0845829027952*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.566058029713*x + -0.0861485259205*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.519391036762*x + -0.0928410129815*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7787986077295 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.62193643511075 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8767732011935 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6956057764083 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8345141816675 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7180845184296 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.86273816015 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.74139595430445 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.52541431265445 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.5358349783497 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.669362791986 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5918492309136 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.6392044787986 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.77537001517098 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.72776177494079 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.619771669669365*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.47444015687955*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.57339046203475 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.632113243580105 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5116087577352 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.58454384064525 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.54805535204115 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.59925362058945 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.637456772545055 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.652267785796*x + -0.09277436198435*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.52466172436537*x + -0.09396923974981*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.6919373380965*x + -0.09929452809875*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.476239852572955 + 0*x + -0.0524981111393145*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.474608149227056 + 0.365177167287945*x + -0.06875109136559*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.058572373596785*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.5098274702395*x + -0.080415748792*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.55857554994815 + 0.3610767696821*x + -0.061534454951435*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.534628832713365*x + -0.0806373364132*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.37194594380655*x + -0.06274479154785*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.558341892320055 + 0.4323289046062*x + -0.0734595934743655*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3721462653705*x + -0.05957674208923*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.37089455401725*x + -0.06232688823925*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.56082547629*x + -0.0841741930285*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7185264079282 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(1.029765876989 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.92602509709689 + -0.773963997471*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.75977711325065*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.69105570573505*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.6909902289927*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.8245713895875*x + 0.083829166254727*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.9496242148358*x + 0.092921151893225*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.712180158669815 + -0.9126669783915*x + 0.092660253843125*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -1.16943191436*x + 0.125110893556165*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -1.2120039088*x + 0.12751219489727*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.71537188626655 + -1.2094400399*x + 0.12997866751755*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.49650887142987*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-2.1016904474 + 0.644322039805*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-2.0863835019 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-2.15313895785 + 0.751962763824*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.56660289276379 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6298852649594 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.56058883024645 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.93029295492*x + -0.08848251822093*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.93304685226*x + -0.082124011054365*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 1.31300288995*x + -0.11003611868269*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.52422133025783*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.57732641850455*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.641339813777815*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.73992341814313*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.6817494638708*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.6161131162948*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.6010977581895*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.39409080031798*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3755516577795*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4109655706004*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.54292430850915*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.60275199457195 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.56904865194735 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.769031518366015 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.71755262540125 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.67584463058336 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.754031931671 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.90098410469678 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8963583325625 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.67942814160285 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5239696937945*x + -0.0663240974723*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.40644109411975*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5696569277365 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.41714119428085*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.59070749323265 + 0.4948267821194*x + -0.063995186227825*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.65988159050645 + 0*x + -0.06470826575195*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.62461722487393 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.664894475684975 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.444719745206362*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.37221068102603*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6047470950353 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.7484065047074*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6145263025442 + 0.6135001051598*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.49381641077076*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.51046611514015*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09476219)+(0.01686738)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #plot mani 2 and 4 are staggered to prevent overlap
  #mani1
  stat_function(fun=function(x){(-0.19162650 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.22 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.24 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.25 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.19162650 + (0.0733269+0.56860900)*x + (-0.00588638-0.05748940)*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#EC1804')

# print(evennessPlot) #export at 1200x1000


#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000


###summary stats from bayesian output --------------------------------------------------------
#gather summary stats needed and relabel them
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

# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_03132017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_03132017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)




###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
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
  scale_y_continuous(limits=c(-0.3, 0.18), breaks=seq(-0.2, 0.2, 0.2)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
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
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Richness Change') +
  annotate('text', x=3.45, y=-0.3, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.4, 0.25), breaks=seq(-0.3, 0.3, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
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





# ###by resource mani - model output--------------------------------------------------------
# #still need to calculate the proportion of chains where x resource response was greater than y resource response
# trtDetail <- expRaw%>%
#   select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
#   group_by(site_code, project_name, community_type, treatment)%>%
#   summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
#   mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))
# 
# chainsTrt <- chainsEquations%>%
#   select(variable, site_code, project_name, community_type, treatment, intercept, linear, quadratic, yr9, yr_final, plot_mani, rrich, anpp, MAT, MAP, min_year, experiment_length, alt_length)%>%
#   left_join(trtDetail)
# 
# resourceMani <- chainsTrt%>%
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
#                                                                                              ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr', 'other'))))))))))))
# 
# #plot model output by resource manipulated at final year of each experiment (varies by experiment) ---------------------------
# meanResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='mean'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.20), name='Mean Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(0, 0.6), xlim=c(1,4)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.60, label='(a)', size=12, hjust='left') +
#   annotate('text', x=1, y=0.42, label='a*', size=10) +
#   annotate('text', x=2, y=0.47, label='ab*', size=10) +
#   annotate('text', x=3, y=0.35, label='b*', size=10) +
#   annotate('text', x=4, y=0.35, label='b*', size=10)
# 
# dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='dispersion'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.01), name='Dispersion Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.05, 0.03), xlim=c(1,4)) +
#   xlab('') +
#   annotate('text', x=0.5, y=0.03, label='(b)', size=12, hjust='left')
# 
# richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='richness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Richness Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.15, 0.06), xlim=c(1,4)) +
#   xlab('') +
#   annotate('text', x=0.5, y=0.06, label='(c)', size=12, hjust='left') +
#   annotate('text', x=1, y=-0.13, label='ab*', size=10) +
#   annotate('text', x=2, y=-0.10, label='a*', size=10) +
#   annotate('text', x=3, y=-0.05, label='b', size=10) +
#   annotate('text', x=4, y=-0.12, label='ab*', size=10)
# 
# evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='evenness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Evenness Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.03, 0.07), xlim=c(1,4)) +
#   xlab('') +
#   annotate('text', x=0.5, y=0.07, label='(d)', size=12, hjust='left') +
#   annotate('text', x=1, y=0.055, label='a*', size=10) +
#   annotate('text', x=3, y=0.05, label='ab*', size=10) +
#   annotate('text', x=4, y=0.02, label='b', size=10)
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600
# 
# 
# ###by magnitude of resource manipulated---------------------------------
# #N addition
# meanNPlotFinal <- ggplot(data=subset(resourceMani, variable=='mean'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=1, label='(a)', size=12, hjust='left')
# 
# dispersionNPlotFinal <- ggplot(data=subset(resourceMani, variable=='dispersion'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=0.5, label='(b)', size=12, hjust='left')
# 
# richnessNPlotFinal <- ggplot(data=subset(resourceMani, variable=='richness'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=2, label='(c)', size=12, hjust='left')
# 
# evennessNPlotFinal <- ggplot(data=subset(resourceMani, variable=='evenness'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(-0.35, 0.5, 0.25), name='Evenness Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=0.5, label='(d)', size=12, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600
# 
# #H2O change
# meanPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='mean'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=1, label='(a)', size=12, hjust='left')
# 
# dispersionPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='dispersion'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=0.5, label='(b)', size=12, hjust='left')
# 
# richnessPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='richness'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=2, label='(c)', size=12, hjust='left')
# 
# evennessPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='evenness'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(-0.35, 0.5, 0.25), name='Evenness Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=0.5, label='(d)', size=12, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600


###by resource mani - raw data--------------------------------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

rawTrt <- rawData%>%
  filter(treatment_year==8|treatment_year==experiment_length)%>%
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
                                                                                             ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr','other'))))))))))))%>%
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
  mutate(n_only=ifelse(n>0&plot_mani==1&site_code!='CEH', 'n only', 'not'))%>%
  mutate(water_only=ifelse(precip!=0&plot_mani==1, 'water only', 'not'))%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type))

#plot raw data by resource manipulated at final year of each experiment (varies by experiment) ---------------------------
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='mean_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.10), name='Mean Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(0, 0.45), xlim=c(1,6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.45, label='(a)', size=12, hjust='left') +
  annotate('text', x=1, y=0.35, label='b*', size=10) +
  annotate('text', x=2, y=0.42, label='a*', size=10) +
  annotate('text', x=3, y=0.335, label='b*', size=10) +
  annotate('text', x=4, y=0.35, label='b*', size=10) +
  annotate('text', x=5, y=0.355, label='b*', size=10) +
  annotate('text', x=6, y=0.33, label='b*', size=10)

dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='dispersion_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Dispersion Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.07, 0.03), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.03, label='(b)', size=12, hjust='left')

richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='S_PC', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.1), name='Richness Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.2, 0.15), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.15, label='(c)', size=12, hjust='left') +
  annotate('text', x=1, y=0.07, label='a', size=10) +
  annotate('text', x=2, y=-0.185, label='b*', size=10) +
  annotate('text', x=3, y=0.14, label='ab', size=10) +
  annotate('text', x=4, y=0.125, label='a*', size=10) +
  annotate('text', x=5, y=-0.115, label='ab', size=10) +
  annotate('text', x=6, y=0.09, label='a', size=10)

evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='SimpEven_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Evenness Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.09, 0.14), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.14, label='(d)', size=12, hjust='left') +
  annotate('text', x=1, y=0.035, label='ab', size=10) +
  annotate('text', x=2, y=0.058, label='a*', size=10) +
  annotate('text', x=3, y=0.128, label='ab', size=10) +
  annotate('text', x=4, y=-0.05, label='b', size=10) +
  annotate('text', x=5, y=-0.055, label='ab', size=10) +
  annotate('text', x=6, y=-0.055, label='b*', size=10)

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 2400 x 1600


###by magnitude of resource manipulated---------------------------------
#N addition
#model
summary(meanModel<-lme(mean_change~n, random=~1|site_proj_comm,data=subset(rawTrt, n_only=='n only')))
summary(dispModel<-lme(dispersion_change~n, random=~1|site_proj_comm,data=subset(rawTrt, n_only=='n only')))
summary(richModel<-lme(S_PC~n, random=~1|site_proj_comm,data=subset(rawTrt, n_only=='n only')))
summary(evenModel<-lme(SimpEven_change~n, random=~1|site_proj_comm,data=subset(rawTrt, n_only=='n only')))

#figs
meanNPlotFinal <- ggplot(data=subset(rawTrt, n_only=='n only'), aes(x=n, y=mean_change)) +
  geom_point(size=5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=0.4, y=1, label='(a)', size=12, hjust='left') +
  geom_segment(aes(x = 0.432, y = 0.2217095, xend = 152.57, yend = 0.6995399), size=3) +
  theme(legend.position="none")
#add text that intercept: t32=7.94, slope: t22=4.15

dispersionNPlotFinal <- ggplot(data=subset(rawTrt, n_only=='n only'), aes(x=n, y=dispersion_change)) +
  geom_point(size=5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=0.4, y=0.4, label='(b)', size=12, hjust='left')
#add text that intercept: n.s., slope: n.s.

richnessNPlotFinal <- ggplot(data=subset(rawTrt, n_only=='n only'), aes(x=n, y=S_PC)) +
  geom_point(size=5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=1.5, label='(c)', size=12, hjust='left') +
  geom_segment(aes(x = 0.432, y = -0.001551951, xend = 152.57, yend = -0.5481047), size=3) +
  theme(legend.position="none")
#add text that intercept: n.s., slope: t22=-2.93

evennessNPlotFinal <- ggplot(data=subset(rawTrt, n_only=='n only'), aes(x=n, y=SimpEven_change)) +
  geom_point(size=5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-0.2, 0.6, 0.2), name='Evenness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.6, label='(d)', size=12, hjust='left') +
  geom_segment(aes(x = 0.432, y = 0.0008551233, xend = 152.57, yend = 0.302005), size=3) +
  theme(legend.position="none")
#add text that intercept: n.s., slope: t22=3.32

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600

#H2O change
#model
summary(meanModelPrecip<-lme(mean_change~precip, random=~1|site_proj_comm,data=subset(rawTrt, water_only=='water only')))
summary(dispModelPrecip<-lme(dispersion_change~precip, random=~1|site_proj_comm,data=subset(rawTrt, water_only=='water only')))
summary(richModelPrecip<-lme(S_PC~precip, random=~1|site_proj_comm,data=subset(rawTrt, water_only=='water only')))
summary(evenModelPrecip<-lme(SimpEven_change~precip, random=~1|site_proj_comm,data=subset(rawTrt, water_only=='water only')))

#figs
meanPrecipPlotFinal <- ggplot(data=subset(rawTrt, water_only=='water only'), aes(x=precip, y=mean_change)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=-80, y=0.65, label='(a)', size=12, hjust='left') +
  geom_segment(aes(x = -80, y = 0.26261887, xend = 116, yend = 0.26261887), size=3) +
  theme(legend.position="none")
#add text that intercept: t24=13.60, slope: n.s.

dispersionPrecipPlotFinal <- ggplot(data=subset(rawTrt, water_only=='water only'), aes(x=precip, y=dispersion_change)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=-80, y=0.2, label='(b)', size=12, hjust='left')
#add text that intercept: n.s., slope: n.s.

richnessPrecipPlotFinal <- ggplot(data=subset(rawTrt, water_only=='water only'), aes(x=precip, y=S_PC)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(-1,2,0.25), name='Richness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.6, label='(c)', size=12, hjust='left')
#add text that intercept: n.s., slope: n.s.

evennessPrecipPlotFinal <- ggplot(data=subset(rawTrt, water_only=='water only'), aes(x=precip, y=SimpEven_change)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(-0.4, 0.5, 0.1), name='Evenness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.2, label='(d)', size=12, hjust='left')
#add text that intercept: n.s., slope: n.s.

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600


# ###for factorial experiments, are responses (agnostic to factor identity) additive?
# #raw data--------------------------------------------------------
# rawFactor <- expRaw%>%
#   select(site_code, project_name, community_type, treatment, factorial)%>%
#   filter(factorial==1)%>%
#   group_by(site_code, project_name, community_type, treatment)%>%
#   unique()%>%
#   left_join(rawData)%>%
#   filter(!is.na(mean_change))%>%
#   filter(treatment_year==9|treatment_year==experiment_length)%>%
#   filter(plot_mani<6)
# 
# #plot factorial responses for each metric
# meanFactorial <- ggplot(barGraphStats(data=rawFactor, variable="mean_change", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity", color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
#   scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
#   xlab('') +
#   annotate('text', x=0.5, y=0.65, label='(a)', size=12, hjust='left')
# 
# dispersionFactorial <- ggplot(barGraphStats(data=rawFactor, variable="dispersion_change", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity", color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
#   scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
#   xlab('') +
#   annotate('text', x=0.5, y=0.38, label='(b)', size=12, hjust='left')
# 
# richnessFactorial <- ggplot(barGraphStats(data=rawFactor, variable="S_PC", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity", color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
#   scale_y_continuous(breaks=seq(-1,2,0.25), name='Richness Change') +
#   xlab('Number of Factors') +
#   annotate('text', x=0.5, y=0.6, label='(c)', size=12, hjust='left')
# 
# evennessFactorial <- ggplot(barGraphStats(data=rawFactor, variable="SimpEven_change", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity", color='black', fill='white') +
#   geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
#   scale_y_continuous(breaks=seq(-0.4, 0.5, 0.1), name='Evenness Change') +
#   xlab('Number of Factors') +
#   annotate('text', x=0.5, y=0.3, label='(d)', size=12, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanFactorial, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionFactorial, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(evennessFactorial, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# print(richnessFactorial, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# #export at 1800 x 1600




#look at number replicates for dispersion results (all factors actually) -- doesn't make a difference --------------------------------------
reps <- read.csv('SpeciesRelativeAbundance_May2017.csv')%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year, plot_id)%>%
  summarise(mean=mean(relcov))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year)%>%
  summarise(rep_num=n())%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(rep_num=mean(rep_num))

dispersionReps <- rawTrt%>%
  left_join(reps, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type))

#model
summary(meanModelDisp<-lme(mean_change~rep_num, random=~1|site_proj_comm,data=dispersionReps))
summary(dispModelDisp<-lme(dispersion_change~rep_num, random=~1|site_proj_comm,data=dispersionReps))
summary(richModelDisp<-lme(S_PC~rep_num, random=~1|site_proj_comm,data=dispersionReps))
summary(evenModelDisp<-lme(SimpEven_change~rep_num, random=~1|site_proj_comm,data=dispersionReps))

#figs
meanRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=mean_change)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Mean Change') +
  scale_x_continuous(breaks=seq(0,20,5)) +
  coord_cartesian(xlim=c(2,20)) +
  annotate('text', x=0, y=1, label='(a)', size=12, hjust='left') +
  geom_segment(aes(x = 3, y = 0.3213704, xend = 18, yend = 0.1860824), size=3) +
  theme(legend.position="none")
#add text that intercept: t210=11.79., slope: t210=-2.16
dispersionRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=dispersion_change)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Dispersion Change') +
  scale_x_continuous(breaks=seq(0,20,5)) +
  coord_cartesian(xlim=c(2,20)) +
  annotate('text', x=0, y=0.35, label='(b)', size=12, hjust='left')
richnessRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=S_PC)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Richness Change') +
  scale_x_continuous(breaks=seq(0,20,5)) +
  coord_cartesian(xlim=c(2,20)) +
  annotate('text', x=0, y=1.4, label='(c)', size=12, hjust='left')
evennessRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=SimpEven_change)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Evenness Change') +
  scale_x_continuous(breaks=seq(0,20,5)) +
  coord_cartesian(xlim=c(2,20)) +
  annotate('text', x=0, y=0.6, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanRepPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionRepPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessRepPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessRepPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600







###look at five factor manipulations for mean change --------------------------------------------------------
#compare any four factor without N to five factor with N - using model data final year--------------------------
expRawMean <- expRaw%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarise(n=mean(n), herb_removal=mean(herb_removal), plant_mani=mean(plant_mani))

meanCompare <- subset(chainsEquations, project_name=='e001'|project_name=='e002'|site_code=='NIN'|site_code=='TRA')%>%
  select(variable, site_code, project_name, community_type, treatment, plot_mani, yr9)%>%
  left_join(expRawMean, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))

#plot without N at four factors, with N at five factors
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&variable=='mean'), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none')
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&herb_removal>0&variable=='mean'), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none')
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plant_mani>0&variable=='mean'), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(labels=c('4 factor\n+manip.', '5 factor\n+manip.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(c) Plant Manipulation Comparison', size=10, hjust='left') +
  theme(legend.position='none')

pushViewport(viewport(layout=grid.layout(1,3)))
print(compareNPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(comparePlantPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(compareHerbPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2400x1200

#compare any four factor without N to five factor with N - using raw data final year--------------------------
expRawMean <- expRaw%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarise(n=mean(n), herb_removal=mean(herb_removal), plant_mani=mean(plant_mani))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, plot_mani, n, herb_removal, plant_mani)

meanCompare <- subset(rawTrt, project_name=='e001'|project_name=='e002'|site_code=='NIN'|site_code=='TRA')%>%
  left_join(expRawMean, all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))


#plot without N at four factors, with N at five factors--------------------------
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6), variable='mean_change', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.4, label='a*', size=10) +
  annotate('text', x=2, y=0.57, label='a*', size=10) +
  annotate('text', x=3, y=0.8, label='b*', size=10)
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6), variable='mean_change', byFactorNames=c('plot_mani', 'herb_removal')), aes(x=interaction(plot_mani, herb_removal), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(limits=c('4.0', '4.1', '5.1'), labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.39, label='*', size=10) +
  annotate('text', x=2, y=0.6, label='*', size=10) +
  annotate('text', x=3, y=0.59, label='*', size=10)
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6), variable='mean_change', byFactorNames=c('plot_mani', 'plant_mani')), aes(x=interaction(plot_mani, plant_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(limits=c('4.0', '4.1', '5.1'), labels=c('4 factor\n-manip', '4 factor\n+manip.', '5 factor\n+manip.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(c) Plant Manipulation Comparison', size=10, hjust='left') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.39, label='*', size=10) +
  annotate('text', x=2, y=0.61, label='*', size=10) +
  annotate('text', x=3, y=0.69, label='*', size=10)

pushViewport(viewport(layout=grid.layout(1,3)))
print(compareNPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(comparePlantPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(compareHerbPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2400x1200


###look at trts with multiple resources --------------------------------------------------------
trtType <- rawTrt%>%
  mutate(trt_type=ifelse(plot_mani==2&n>0&drought<0, 'N+drought', ifelse(plot_mani==2&n>0&irrigation>0, 'N+irr', ifelse(plot_mani==2&n>0&p>0, 'N+P', ifelse(plot_mani==2&p>0&k>0, 'P+K', ifelse(plot_mani==2&n>0&CO2>0, 'N+CO2', ifelse(plot_mani==2&CO2>0&irrigation>0, 'CO2+irr', ifelse(plot_mani==3&n>0&p>0&k>0, 'N+P+K', ifelse(plot_mani==3&n>0&CO2>0&irrigation>0, 'N+CO2+irr', ifelse(plot_mani==4&n>0&p>0&k>0&irrigation>0, 'N+P+K+irr', ifelse(plot_mani==1&n>0, 'N', ifelse(plot_mani==1&p>0, 'P', ifelse(plot_mani==1&irrigation>0, 'irr', ifelse(plot_mani==1&drought<0, 'drought', ifelse(plot_mani==1&CO2>0, 'CO2', 'other')))))))))))))))%>%
  mutate(disp_abs=abs(dispersion_change), rich_abs=abs(S_PC), even_abs=abs(SimpEven_change))

meanChangeInteractions <- ggplot(barGraphStats(data=subset(trtType, trt_type!='other'), variable="mean_change", byFactorNames=c("trt_type")), aes(x=trt_type, y=mean))+
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_x_discrete(limits=c('N', 'P', 'drought', 'irr', 'CO2', 'N+drought', 'N+irr', 'N+P', 'P+K', 'N+CO2', 'CO2+irr', 'N+P+K', 'N+CO2+irr', 'N+P+K+irr')) +
  theme(axis.text.x  = element_blank()) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=11.5) +
  geom_vline(xintercept=13.5) +
  ylab('Mean Change')

dispersionChangeInteractions <- ggplot(barGraphStats(data=subset(trtType, trt_type!='other'), variable="disp_abs", byFactorNames=c("trt_type")), aes(x=trt_type, y=mean))+
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_x_discrete(limits=c('N', 'P', 'drought', 'irr', 'CO2', 'N+drought', 'N+irr', 'N+P', 'P+K', 'N+CO2', 'CO2+irr', 'N+P+K', 'N+CO2+irr', 'N+P+K+irr')) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  theme(axis.text.x  = element_blank()) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=11.5) +
  geom_vline(xintercept=13.5) +
  ylab('Abs Dispersion Change')

richnessChangeInteractions <- ggplot(barGraphStats(data=subset(trtType, trt_type!='other'), variable="rich_abs", byFactorNames=c("trt_type")), aes(x=trt_type, y=mean))+
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_x_discrete(limits=c('N', 'P', 'drought', 'irr', 'CO2', 'N+drought', 'N+irr', 'N+P', 'P+K', 'N+CO2', 'CO2+irr', 'N+P+K', 'N+CO2+irr', 'N+P+K+irr')) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  theme(axis.text.x  = element_blank()) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=11.5) +
  geom_vline(xintercept=13.5) +
  ylab('Abs Richness Change')

evennessChangeInteractions <- ggplot(barGraphStats(data=subset(trtType, trt_type!='other'), variable="even_abs", byFactorNames=c("trt_type")), aes(x=trt_type, y=mean))+
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.2) +
  scale_x_discrete(limits=c('N', 'P', 'drought', 'irr', 'CO2', 'N+drought', 'N+irr', 'N+P', 'P+K', 'N+CO2', 'CO2+irr', 'N+P+K', 'N+CO2+irr', 'N+P+K+irr')) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  theme(axis.text.x  = element_blank()) +
  geom_vline(xintercept=5.5) +
  geom_vline(xintercept=11.5) +
  geom_vline(xintercept=13.5) +
  ylab('Abs Evenness Change')


pushViewport(viewport(layout=grid.layout(4,1)))
print(meanChangeInteractions, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionChangeInteractions, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(richnessChangeInteractions, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(evennessChangeInteractions, vp=viewport(layout.pos.row = 4, layout.pos.col = 1))
#export at 1800 x 1600
