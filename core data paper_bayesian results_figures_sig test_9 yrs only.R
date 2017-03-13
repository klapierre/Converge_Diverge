library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

setwd("C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

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
expRaw <- read.csv('ExperimentInformation_Dec2016.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('ForBayesianAnalysis_9yr_Dec2016.csv')

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

#raw chains data --------------------------------------------------------
chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_cholesky_3.csv', comment.char='#')
chains4 <- chains4[-1:-5000,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)


#density plot of chains --------------------------------------------------------
plot(density(chainsCommunity$mu.1.1))
plot(density(chainsCommunity$mu.1.2))
plot(density(chainsCommunity$mu.1.3))


#get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, 
         #plot_mani intercepts (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1.1.1, U.2.1.1, U.3.1.1, U.4.1.1,
         U.1.2.1, U.2.2.1, U.3.2.1, U.4.2.1,
         U.1.3.1, U.2.3.1, U.3.3.1, U.4.3.1,
         U.1.4.1, U.2.4.1, U.3.4.1, U.4.4.1,
         #plot_mani linear slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1.1.2, U.2.1.2, U.3.1.2, U.4.1.2,
         U.1.2.2, U.2.2.2, U.3.2.2, U.4.2.2,
         U.1.3.2, U.2.3.2, U.3.3.2, U.4.3.2,
         U.1.4.2, U.2.4.2, U.3.4.2, U.4.4.2,
         #plot_mani quad slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
         U.1.1.3, U.2.1.3, U.3.1.3, U.4.1.3,
         U.1.2.3, U.2.2.3, U.3.2.3, U.4.2.3,
         U.1.3.3, U.2.3.3, U.3.3.3, U.4.3.3,
         U.1.4.3, U.2.4.3, U.3.4.3, U.4.4.3,
         #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
         D.1.1.1, D.2.1.1, D.3.1.1, D.4.1.1,
         D.1.1.2, D.2.1.2, D.3.1.2, D.4.1.2,
         D.1.1.3, D.2.1.3, D.3.1.3, D.4.1.3,
         #richness intercept, linear, and quad slopes (center digit): 2=richness
         D.1.2.1, D.2.2.1, D.3.2.1, D.4.2.1,
         D.1.2.2, D.2.2.2, D.3.2.2, D.4.2.2,
         D.1.2.3, D.2.2.3, D.3.2.3, D.4.2.3,
         #MAP intercept, linear, and quad slopes (center digit): 1=MAP
         E.1.1.1, E.2.1.1, E.3.1.1, E.4.1.1,
         E.1.1.2, E.2.1.2, E.3.1.2, E.4.1.2,
         E.1.1.3, E.2.1.3, E.3.1.3, E.4.1.3,
         #MAT intercept, linear, and quad slopes (center digit): 2=MAT
         E.1.2.1, E.2.2.1, E.3.2.1, E.4.2.1,
         E.1.2.2, E.2.2.2, E.3.2.2, E.4.2.2,
         E.1.2.3, E.2.2.3, E.3.2.3, E.4.2.3,
         #overall intercept, linear, and quad slopes
         mu.1.1, mu.2.1, mu.3.1, mu.4.1,
         mu.1.2, mu.2.2, mu.3.2, mu.4.2,
         mu.1.3, mu.2.3, mu.3.3, mu.4.3)%>%
  gather(key=parameter, value=value, U.1.1.1:mu.4.3)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))



#gather the intercepts, linear slopes, and quadratic slopes for all treatments --------------------------------------------------------
#numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
#variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
#parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
#set any that are not significant (CI overlaps 0) as 0

#get mean parameter values across all runs for each experiment, treatment, etc
chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,12172:17235]))%>% #may need to delete original four chains dataframes to get this to work
  add_rownames('parameter')
names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 12172:17235])'] <- 'mean'
#get sd of parameter values across all runs for each experiment, treatment, etc
chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,12172:17235]))
names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 12172:17235])'] <- 'sd'

chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
  #split names into parts
  separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
  select(-B)%>%
  #rename parts to be more clear
  mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
         parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
         id=as.integer(id))%>%
  #if 95% confidence interval overlaps 0, then set mean to 0
  mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
  #spread by variable
  select(variable, id, parameter, mean)%>%
  spread(key=parameter, value=mean)

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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.657442823605295 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.55771095794565 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.790312229892 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.057372042471455*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.63469624632414 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6708266664662 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.80042527721565 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.90554047621985 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.1965505809 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.8318747363 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4136183434196*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7579641981388 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4343147872223*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.40902394215 + -0.42341762914575*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.887163073428 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04709100823598*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.337810538111*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7779624379276 + 0.38833324960852*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5202541935285*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.533996787465*x + -0.06161208345553*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8735192727382 + -0.47267243144815*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6021760213715*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.488275445664*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.7372195288563 + -0.43209328441995*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.992081653365 + -0.33722307894663*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.940065754505 + -0.4363325293945*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.10318979694 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.8041759221039 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.75298531871555 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.68042609026815 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.62107787044655 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.40567163515315*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.59137681585715 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.85070661083255 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.68376517395755 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0345251319262585*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03981295669433*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0349243767412157*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.23740433568*x + -0.038594125619558*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.58293819045 + 0.92210810455*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.2142653915 + 0.858244469475*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.6832231968 + 0.850817143415*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.5640353523 + 0.937181674225*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.265701047 + 0.7872133299455*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.36612152015 + 1.031838637555*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.90994949735 + 0.85892738673*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.23974779227365 + 1.0089132254*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.7197172393 + 0.933873780995*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.2345271825 + 0.81199587453*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.0751416159145 + 1.0907217838*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7822986216345 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.37025265738045*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7585390664745 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.80917757076985 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.78601541759 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8313921433613 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.89790809260755 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.67746516394025 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.69772419369473 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.67988073978975 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8128797907106 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7958562383314 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.9844944493542 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.737156883401385 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.657114007982445 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.156525203345 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.047260101513555*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.163439440115 + 0.6142167444*x + -0.082494122195*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.072089077535 + 0.4280360217892*x + -0.0805663774728*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.820706838456 + 0.38282957116115*x + -0.066844929449*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.01279516261 + 0.458766505376855*x + -0.066941204693*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.12326691519 + 0.543086471395*x + -0.08737341455*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.146920168495 + 0.455491110593355*x + -0.080672881357*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.06331524727*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04177531892945*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.038245393985935*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7157094870954 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8559718681595 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4401516962308*x + 0.091394925749*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4849498888963*x + 0.06623948591345*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3549564680971*x + 0.088226507695*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.60672827855905 + -0.62556353197*x + 0.08658100577*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.018417545775 + -0.590253739665*x + 0.067107181599*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.74323738145*x + 0.09006407014*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.183455726 + -0.455160336176*x + 0.072881286544*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.60398075551811 + -0.65800409635*x + 0.082441887535*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.756552368814 + -0.5629766218005*x + 0.0808602113635*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.70758411663855 + -0.596358710135*x + 0.095207328965*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.85006998698035 + -0.4965579965353*x + 0.07547826595775*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.63883119344*x + 0.086787250085*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.71415300368015 + -0.60893174465*x + 0.081431784455*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8084495969858 + -0.586547398565*x + 0.072393596484*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.70217208444837 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.86437705383 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.890922390236 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.24233667856 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.1649062954655 + 0.44634407171995*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.60581688079215 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.31089597123 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.9738634289595 + 0.50281154581015*x + -0.06429442829854*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0.33282335726005*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.74631846260965 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.85272397203635 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6960531132731 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.9157323461205 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.90554888444425 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8525188183778 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.105691220955 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.5182031507945*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.984289159223 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38065946306425*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.47946453625295*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3520073804151*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(1.0162404553 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.80213836502455 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.86701522343415 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.0914814155168 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6704376339514 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-1.25869828105 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.893045544099 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9017718791605 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.76984509224285 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8009023795355 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.57364814025 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.3495770064 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.862931061464 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.13472116055 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5834755228865 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.42948208286475*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7092779445521 + 0*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4743388065892*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.361905479375355*x + 0*x^2)*(0.08790049)+(-0.000407613)}, size=0.5, xlim=c(0,2), colour='grey') +
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
  ylab('Proportion Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.052500725680245*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.053720770629994*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.604669627567 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7119874327321 + 0*x + 0.053225233406455*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.61889914793265 + 0*x + 0.052964003456365*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.9382452189805 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.058355689471764*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.52178273273539 + 0*x + 0.063457684700065*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6709793516807 + 0*x + 0.05297467529186*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.60844813458901 + 0*x + 0.056062412380698*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5823934732531 + 0*x + 0.058425181299365*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.731933671320915 + 0*x + 0.061051417310735*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.545686114581515 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.600765941543 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.66065652420918 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.72111158436145 + -0.3365595382697*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3181202123263*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4370084022505*x + 0.0370906806338295*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.684306782521005 + -0.291506713236955*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.0955339924335 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.69134068838205 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.96933282003755 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.069962301335 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.75382636638075 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5670644575886*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.77810297561115 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(1.065512892008 + -0.65970009795*x + 0.051318374299*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.6419003162381 + -1.25320598305*x + 0.100301814165*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.3116676041*x + 0.113341521725*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.8682660647*x + 0.06787613102317*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.5294203788 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(1.5376194033 + -0.9650724209*x + 0.0671937742805*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.96263139905*x + 0.088841052375*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.98318783645*x + 0.088377491685*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.5520064455735*x + 0.0586241736675*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.93226562695*x + 0.098945142865*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.7148623819355 + -0.42954241264295*x + 0.05590312604182*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.59429128138596 + -0.474051136431*x + 0.04322824965388*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.25178181765*x + 0.113049830995*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.13025070785*x + 0.09483104751*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.608591622169385 + -0.3760658572105*x + 0.046343895274635*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.0493483781*x + 0.097194085535*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.24564072235*x + 0.10922013392*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.33064247090569*x + 0.044360744309675*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.89651104325*x + 0.08999293198*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -1.0230700936*x + 0.10102485558*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.6034714692813 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6891956637691 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.222563331035 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.71606182122395 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.47675553663905*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(1.2053201958 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.99456635943 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.9126079217155 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7493208195502 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.5431971584 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.744961088602 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.36477356805 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.95381765999 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.48127920005 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.695494465489 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.3863347034 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.60757473625235 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.159097428215 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4923964272897*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.50705523642595*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0395436969298715*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.473562438108*x + -0.06726917855295*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3485608907598*x + -0.0650587881608*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.055278374697395*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.04924721484689*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(1.09107637695 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.9434157261825 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.44205548367174*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.36079220239165*x + -0.043038905473715*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.33742501081415*x + -0.0455854749477*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.289364347055*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.29625180885077*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.274085383489665*x + -0.039902670970815*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.326177237622545*x + -0.0494091040823*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.257721993681246*x + -0.031669633834555*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.30068052858001*x + -0.0369109297830615*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.03170442490302*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.281829683257415*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3538936318484*x + -0.043114892596445*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.24575813372062*x + -0.0358756559443725*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.03020412730093*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.781212935955*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.7583844937319*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.68565916020465*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.76484403442155*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.69466048068165*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.64368342677695*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4632579230268*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.29461644333779*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3415846246601*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.60934342014823 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 1.31322693094015*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 1.2318801074065*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7763664299612 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.52705523507805*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.58937325799635*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5286150654389*x + -0.052672873888455*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8749867156714 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.5654271666505*x + -0.05852840985879*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.856523746876 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7783262941175 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.8226685305075 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.822005318538856 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.31948518890595*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5197281236405 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.38664181636935*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.36464874499113*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0384706116446215*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.037534651118735*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.575012182214707 + 0*x + -0.04187530732092*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(1.0606560420055 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.22740690754 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.138034600355 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.229078646155 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.996977617484 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.9426322806765 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.68952217778376 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4205437723941*x + 0.05424925390994*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.700280031606755 + -0.4664737093325*x + 0.06010632655985*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.57965823513925 + 0*x + 0.0476181753023005*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0467137791247*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.247553493774545*x + 0.04553869483361*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.43418397327795*x + 0.064931562886*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.26569759776276*x + 0.03913943104054*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.37031315038985*x + 0.05697863700776*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.523305999141075 + -0.36954012055045*x + 0.064631404626*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0558710423435*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.03741452358436*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.04723849716805*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.04557486806822*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.6142385669792 + -0.52694284964215*x + 0.058444915469595*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8988182791305 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.43146848846645*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.78679700085835 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.6199326485035*x + 0.069272749559825*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.63474081608735*x + 0.071072423922525*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.70760293153205 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.9111281805595 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.822140399050775 + -0.61332232902855*x + 0.0685051832490345*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.8182409502085 + -0.5516850439787*x + 0.06470399765685*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.60003238875 + -0.554336418933795*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.632764641155 + -2.38334280005955*x + 0.86712843115*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(2.32586489325 + 0*x + 0.81166994528771*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(2.30977549435 + -1.80149796316365*x + 0.8163838465643*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.0595943513 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.6509192656773*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.47170656046538*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.38628407991624*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.59317668763465*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.6253704551275*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5572575562286*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.75852615135685 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.770493619459395 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.43898400035105*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.40533545277015*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.32426552806891*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.37831771302815*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(1.0052484282 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.53589640701277 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.563270355128285 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.956199195253 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.60363560721555 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0.68422115271255 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7528920050958 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0.58557235014315 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.60222617728139 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.657240722046265 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.75762403256095 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.78536591751825 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.927804865313 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.65439312887911 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.511542998941675*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.946152104148 + 0.5779647476425*x + -0.07616841231025*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4375012117428*x + -0.059326889065569*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.777834160988755 + 0.517135206617*x + -0.0656142898325*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.9363308807215 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.767912509129 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.968938597367 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.38006061135355*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7255134447973 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.58571206872575 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.72936022353455 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.76351867818525 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.161810812605385*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.7573243263717 + 0*x + 0.159640615054465*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.16755140944097*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.443117945152745*x + -0.11828345047559*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.44271262900045*x + -0.12597805930595*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.70297505349548*x + -0.14946251576105*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.7170046174705 + 0.50557943441689*x + -0.13572609570035*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.899635291505 + 0.4710126578241*x + -0.13734951157645*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.12099889694403*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.12388321785341*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.11766924315015*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.9229009926 + 0.8292686646975*x + -0.17379465041715*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.59793095410065 + 0.7037165466795*x + -0.151819340077*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.114528654305806*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.66320734953785 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.40955049705 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9066701714328 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.23759067875 + 0*x + 0*x^2)*(0.2154907)+(-0.0546117)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #plot mani 1-4 staggered to prevent overlap
  #mani1
  stat_function(fun=function(x){(0.30632550 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(0.32 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(0.34 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(0.36 + 0*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(0.38 + (-0.0670179-0.40687600)*x + 0*x^2)*0.2154907 + -0.0546117}, size=3, xlim=c(0,8), colour='#EC1804')

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
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6262137145461 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.44791070858958*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6051760813038 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.51942709207485 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.609400058802605 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6556295452597 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.50218452545715*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.8122071118208 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.878173560685 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.05780014575 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7728725684575 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.010826571584 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3797075573566*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4869975519824*x + -0.05190305007372*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.135157479376 + 0.6507617902705*x + -0.058273867138875*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.3763511384 + 1.35265025755*x + -0.120495927645*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-1.117954816365 + 1.882056908*x + -0.15884103651*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7820606732003 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.61284291952*x + -0.05731822191215*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.418032407072235*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7252829777496 + 0.376417370607465*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.493008281294365*x + -0.044498920499365*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.9858780515*x + -0.175521407605*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 1.2293092911*x + -0.11279284892*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.72939387641445 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.0730174296485 + 0.83623121974*x + -0.07921693553155*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.8669123140193 + 1.12602281405*x + -0.08642683578275*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.82787638388075 + 0.53601337600005*x + -0.04837418802648*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.8325054212*x + -0.07167478965277*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.47066436153044 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5289472630016 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.47251022175565*x + -0.05491627547123*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.43359565073842*x + -0.0759093067239*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.374760240290945*x + -0.0702240658066205*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4565303069791*x + -0.07890789787725*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.48513373295285*x + -0.07559510553615*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4368375925951*x + -0.0789997217363*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.52897174310245*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.746770463839 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.60445343263695 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8364652039048 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.51947612897529 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.670471525619845 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7973310371505 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6914800752515 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8222107284607 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.71512616352485 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.61335796125755 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.53892746336925 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.66159700958325 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.438114088688587*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.79490150560375 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.56836590636595*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3817962270332*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5494058219869*x + -0.049562128730985*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.38058543914425*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.370945684942755*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.32886526594275*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5690788992522 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.59044707036838 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.51704579723375 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.55491961935685 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5488211091699 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.587266348308725 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.592536340077945 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.406524110216772*x + -0.049604198445765*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.375753894593545*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.37469248092875*x + -0.0478008944469185*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.038840217190425*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.9342808402565 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.84521780625335 + -0.6961499796335*x + 0.06567857831117*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.69017931334195*x + 0.071975336098875*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.62354525732389*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.590443482599415*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.6771718181454*x + 0.069758884632685*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.8122939868106*x + 0.081466074323575*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.7393622377021*x + 0.075485341397121*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.989625007537145*x + 0.1034697051059*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -1.03444529902*x + 0.10716508352885*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.9966434618835*x + 0.1044776015926*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.9953370953 + 0.57284192769835*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.97886323855 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-2.044502377 + 0.689565913695215*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.580836635628105 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7276998328982 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.57133140781615 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.59295539651195*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4569126587435*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8055399853771 + 1.18374435725*x + -0.101117191295865*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.8174004723845*x + -0.071327295421055*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.859954344201*x + -0.0869084168665*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.87592521483*x + -0.08885160691914*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.599545461949*x + -0.06694513362484*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.53265407566348*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.6081527127345*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.6500379885525*x + -0.07336320124268*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.40824234759005*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.34806282406085*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.36134449578965*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.397699879077*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3674652029652*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.51251736494955*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.537955639044795 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7445122613507 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.67728187675525 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.68123964475015 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6419563787603 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7145733718074 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8554071754633 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.85021000146935 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.64903799951995 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.42101055173735*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.523367486945595 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.522183436977 + 0.39089989755325*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.593769500393055 + 0*x + -0.057700970997245*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5941902691266 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6314777443888 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.38443466058754*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.54381557218195 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.54769274147525 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.894545802945*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.73769722283215 + 0.754399415107*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.540607169743475*x + -0.06620800737132*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #plot mani 2 and 4 are staggered to prevent overlap
  #mani1
  stat_function(fun=function(x){(-0.18698100 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.2 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.22 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.24 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.18698100 + (0.04984795+0.71206950)*x + (-0.00588638-0.05748940)*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,8), colour='#EC1804')

# print(evennessPlot) #export at 1200x1000


#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000




###by resource mani --------------------------------------------------------
#still need to calculate the proportion of chains where x resource response was greater than y resource response
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

chainsTrt <- chainsEquations%>%
  select(variable, site_code, project_name, community_type, treatment, intercept, linear, quadratic, yr9, plot_mani, rrich, anpp, MAT, MAP, min_year, experiment_length, alt_length)%>%
  left_join(trtDetail)

resourceMani <- chainsTrt%>%
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
                                                                                             ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr', 'other'))))))))))))

#plot by resource manipulated at year 9 --------------------------------------------------------
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='mean'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.20), name='Mean Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(0, 0.6), xlim=c(1,4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.60, label='(a)', size=12, hjust='left') +
  annotate('text', x=1, y=0.43, label='*', size=10) +
  annotate('text', x=2, y=0.455, label='*', size=10) +
  annotate('text', x=3, y=0.39, label='*', size=10) +
  annotate('text', x=4, y=0.41, label='*', size=10)

dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='dispersion'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Dispersion Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.05, 0.17), xlim=c(1,4)) +
  xlab('') +
  annotate('text', x=0.5, y=0.17, label='(b)', size=12, hjust='left') +
  annotate('text', x=4, y=0.16, label='*', size=10)

richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='richness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.2, 0.28), xlim=c(1,4)) +
  xlab('') +
  annotate('text', x=0.5, y=0.28, label='(c)', size=12, hjust='left') +
  annotate('text', x=2, y=-0.195, label='*', size=10)

evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='evenness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Evenness Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.04, 0.07), xlim=c(1,4)) +
  xlab('') +
  annotate('text', x=0.5, y=0.07, label='(d)', size=12, hjust='left') +
  annotate('text', x=1, y=0.06, label='*', size=10) +
  annotate('text', x=3, y=0.055, label='*', size=10)

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600






###summary stats from bayesian output --------------------------------------------------------
#gather summary stats needed and relabel them
chainsCommunitySummary <- chainsCommunity%>%
  select(#plot_mani intercepts (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
    U.1.1.1, U.2.1.1, U.3.1.1, U.4.1.1,
    U.1.2.1, U.2.2.1, U.3.2.1, U.4.2.1,
    U.1.3.1, U.2.3.1, U.3.3.1, U.4.3.1,
    U.1.4.1, U.2.4.1, U.3.4.1, U.4.4.1,
    #plot_mani linear slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
    U.1.1.2, U.2.1.2, U.3.1.2, U.4.1.2,
    U.1.2.2, U.2.2.2, U.3.2.2, U.4.2.2,
    U.1.3.2, U.2.3.2, U.3.3.2, U.4.3.2,
    U.1.4.2, U.2.4.2, U.3.4.2, U.4.4.2,
    #plot_mani quad slopes (center digit): 1=plot mani 2, 2=plot mani 3, 3=plot mani 4, 4=plot mani 5
    U.1.1.3, U.2.1.3, U.3.1.3, U.4.1.3,
    U.1.2.3, U.2.2.3, U.3.2.3, U.4.2.3,
    U.1.3.3, U.2.3.3, U.3.3.3, U.4.3.3,
    U.1.4.3, U.2.4.3, U.3.4.3, U.4.4.3,
    #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
    D.1.1.1, D.2.1.1, D.3.1.1, D.4.1.1,
    D.1.1.2, D.2.1.2, D.3.1.2, D.4.1.2,
    D.1.1.3, D.2.1.3, D.3.1.3, D.4.1.3,
    #richness intercept, linear, and quad slopes (center digit): 2=richness
    D.1.2.1, D.2.2.1, D.3.2.1, D.4.2.1,
    D.1.2.2, D.2.2.2, D.3.2.2, D.4.2.2,
    D.1.2.3, D.2.2.3, D.3.2.3, D.4.2.3,
    #MAP intercept, linear, and quad slopes (center digit): 1=MAP
    E.1.1.1, E.2.1.1, E.3.1.1, E.4.1.1,
    E.1.1.2, E.2.1.2, E.3.1.2, E.4.1.2,
    E.1.1.3, E.2.1.3, E.3.1.3, E.4.1.3,
    #MAT intercept, linear, and quad slopes (center digit): 2=MAT
    E.1.2.1, E.2.2.1, E.3.2.1, E.4.2.1,
    E.1.2.2, E.2.2.2, E.3.2.2, E.4.2.2,
    E.1.2.3, E.2.2.3, E.3.2.3, E.4.2.3,
    #overall intercept, linear, and quad slopes
    mu.1.1, mu.2.1, mu.3.1, mu.4.1,
    mu.1.2, mu.2.2, mu.3.2, mu.4.2,
    mu.1.3, mu.2.3, mu.3.3, mu.4.3)%>%
  gather(key=parameter, value=value, U.1.1.1:mu.4.3)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(CI=sd*2)%>%
  separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
  mutate(parameter=ifelse(level=='mu', predictor, parameter), predictor=ifelse(level=='mu', 'overall', predictor))%>%
  #rename parts to be more clear
  mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
         parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
         predictor=ifelse(level=='D'&predictor==1, 'ANPP', ifelse(level=='D'&predictor==2, 'rrich', ifelse(level=='E'&predictor==1, 'MAP', ifelse(level=='E'&predictor==2, 'MAT', ifelse(level=='U'&predictor==1, 'plot mani 2', ifelse(level=='U'&predictor==2, 'plot mani 3', ifelse(level=='U'&predictor==3, 'plot mani 4', ifelse(level=='U'&predictor==4, 'plot mani 5', 'overall')))))))))%>%
  select(level, parameter, variable, predictor, predictor, median, sd, CI)

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)


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


#overall responses --------------------------------------------------------
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
  ggtitle('Proportion\nRichness Change') +
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



###look for patterns of spp appearance/disappearance -- no clear patterns, probably because just the few CDR examples that are long term enough to see the pattern--------------------------------------------------------
relAbund <- read.csv('SpeciesRelativeAbundance_Dec2016.csv')%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_trt=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
  #get rid of duplicate species within a plot and year in the dataset; once we contact the dataowners, this step will no longer be needed
  group_by(exp_trt, site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species)%>%
  summarise(relcov=mean(relcov))%>%
  filter(exp_trt!='NIN::herbdiv::0::5F' & site_code!='GVN')

expinfo<-read.csv('ExperimentInformation_Mar2016.csv')%>%
  mutate(exp_trt=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
  select(exp_trt, plot_mani, calendar_year)

relAbundYear<-merge(relAbund, expinfo, by=c("exp_trt","calendar_year"), all=F)

#make a new dataframe with just the label
exp_trt=relAbundYear%>%
  select(exp_trt)%>%
  unique()

#make a new dataframe to collect the turnover metrics
turnoverAll=data.frame(row.names=1)

for(i in 1:length(relAbundYear$exp_trt)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=relAbundYear[relAbundYear$exp_trt==as.character(exp_trt$exp_trt[i]),]%>%
    select(exp_trt, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
    #get just first and last year of study
    filter(calendar_year==min(calendar_year)|calendar_year==max(calendar_year))
  
  #need this to keep track of plot mani
  labels=subset%>%
    select(exp_trt, plot_mani, calendar_year)%>%
    unique()
  
  #calculate disappearance
  disappearance=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', replicate.var=NA, metric='disappearance')%>%
    group_by(calendar_year)%>%
    summarise(disappearance=mean(disappearance))
  
  #calculate appearance
  appearance=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', replicate.var=NA, metric='appearance')%>%
    group_by(calendar_year)%>%
    summarise(appearance=mean(appearance))
  
  #calculate turnover
  total=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', replicate.var=NA, metric='total')%>%
    group_by(calendar_year)%>%
    summarise(turnover=mean(total))
  
  #merging back with labels to get back plot_mani
  turnover=labels%>%
    left_join(disappearance, by='calendar_year')%>%
    left_join(appearance, by='calendar_year')%>%
    left_join(total, by='calendar_year')%>%
    filter(calendar_year==max(calendar_year))%>%
    select(exp_trt, plot_mani, appearance, disappearance, turnover)
  
  #pasting variables into the dataframe made for this analysis
  turnoverAll=rbind(turnover, turnoverAll)
}

turnoverCtl <- turnoverAll%>%
  filter(plot_mani==0)%>%
  separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)%>%
  select(site_code, project_name, community_type, appearance, disappearance, turnover)
names(turnoverCtl)[names(turnoverCtl)=='appearance'] <- 'appearance_ctl'
names(turnoverCtl)[names(turnoverCtl)=='disappearance'] <- 'disappearance_ctl'
names(turnoverCtl)[names(turnoverCtl)=='turnover'] <- 'turnover_ctl'

turnoverDiff <- turnoverAll%>%
  mutate(trt=ifelse(plot_mani==0, 'ctl', 'trt'))%>%
  separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)%>%
  filter(trt!='ctl')%>%
  left_join(turnoverCtl, by=c('site_code', 'project_name', 'community_type'))%>%
  mutate(appearance_diff=appearance-appearance_ctl, disappearance_diff=disappearance-disappearance_ctl, turnover_diff=turnover-turnover_ctl)

plot(turnoverDiff$plot_mani, turnoverDiff$appearance_diff)
plot(turnoverDiff$plot_mani, turnoverDiff$disappearance_diff)
plot(turnoverDiff$plot_mani, turnoverDiff$turnover_diff)
plot(turnoverDiff$plot_mani, turnoverDiff$appearance)
plot(turnoverDiff$plot_mani, turnoverDiff$disappearance)
plot(turnoverDiff$plot_mani, turnoverDiff$turnover)

summary(glm(turnover~as.factor(plot_mani), data=turnoverDiff))
lsmeans(glm(turnover~as.factor(plot_mani), data=turnoverDiff), 'plot_mani')

turnoverPlot <- ggplot(data=turnoverDiff, aes(x=as.factor(plot_mani), y=turnover)) +
  geom_boxplot() +
  xlab('') +
  ylab('Species Turnover')

appearancePlot <- ggplot(data=turnoverDiff, aes(x=as.factor(plot_mani), y=appearance)) +
  geom_boxplot() +
  xlab('') +
  ylab('Species Appearance')

disappearancePlot <- ggplot(data=turnoverDiff, aes(x=as.factor(plot_mani), y=disappearance)) +
  geom_boxplot() +
  xlab('Number of Factors Manipulated') +
  ylab('Species Disappearance')

pushViewport(viewport(layout=grid.layout(3,1)))
print(turnoverPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(appearancePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(disappearancePlot, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
#export at 1400x3000

#export at 1400x1000

# #compares richness change to turnover metrics
# turnoverRichness <- richness4%>%
#   left_join(turnoverDiff, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
#   select(site_code, project_name, community_type, treatment, experiment_length, plot_mani, intercept, slope, quad, min_year, nutrients, water, carbon, precip, alt_length, yr9, yr20, appearance, disappearance, appearance_diff, disappearance_diff)%>%
#   filter(slope<0, quad>0)
# 
# plot(turnoverRichness$quad, turnoverRichness$appearance_diff)
# plot(turnoverRichness$quad, turnoverRichness$disappearance_diff)
# plot(turnoverRichness$quad, turnoverRichness$appearance)
# plot(turnoverRichness$quad, turnoverRichness$disappearance)
# plot(turnoverRichness$yr9, turnoverRichness$appearance_diff)
# plot(turnoverRichness$yr9, turnoverRichness$disappearance_diff)




# ###look at spp comp of five factor manipulations to find patterns of immigration or loss of dominant spp
# relAbundFive <- relAbundYear
# 
# #make a new dataframe with just the label
# expTrtYear=relAbundFive%>%
#   select(exp_trt)%>%
#   unique()
# 
# #make a new dataframe to collect the turnover metrics
# relAbundFiveYear=data.frame(row.names=1)
# 
# for(i in 1:length(expTrtYear$exp_trt)) {
# 
#   #creates a dataset for each unique year, trt, exp combo
#   subset=relAbundFive[relAbundFive$exp_trt==as.character(expTrtYear$exp_trt[i]),]%>%
#     select(exp_trt, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
#     group_by(exp_trt, calendar_year, treatment, plot_mani, genus_species)%>%
#     summarise(relcov=mean(relcov))%>%
#     ungroup()%>%
#     #get just first and last year of study
#     filter(calendar_year==min(calendar_year)|calendar_year==max(calendar_year))%>%
#     mutate(time=ifelse(calendar_year==min(calendar_year), 'first', 'last'))%>%
#     select(-calendar_year, -treatment)%>%
#     spread(key=time, value=relcov, fill=0)
# 
#   #pasting variables into the dataframe made for this analysis
#   relAbundFiveYear=rbind(subset, relAbundFiveYear)
# }
# 
# relAbundFiveYear <- relAbundFiveYear%>%
#   separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)
# 
# relAbundFiveYearRich <- richness4%>%
#   left_join(relAbundFiveYear, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
#   select(site_code, project_name, community_type, treatment, plot_mani, genus_species, first, last, experiment_length, intercept, slope, quad, nutrients, water, carbon, precip, alt_length, yr9)%>%
#   filter(plot_mani>4)
# 
# ggplot(data=relAbundFiveYearRich, aes(x=first, y=last, colour=quad)) +
#   geom_point() +
#   scale_colour_gradientn(colours=rainbow(4))





#look at number replicates for dispersion results -- doesn't make a difference --------------------------------------------------------
reps <- relAbund%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year, plot_id)%>%
  summarise(mean=mean(relcov))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year)%>%
  summarise(rep_num=n())%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(rep_num=mean(rep_num))

dispersionReps <- subset(resourceMani, variable=='dispersion')%>%
  left_join(reps, by=c('site_code', 'project_name', 'community_type', 'treatment'), all=F)%>%
  #filter out lovegrass mistake trt code until fixed in main dataset
  filter(rep_num>1)

# ggplot(data=dispersionReps, aes(x=rep_num, y=intercept, colour=plot_mani)) +
#   geom_point()
# ggplot(data=dispersionReps, aes(x=rep_num, y=slope, colour=plot_mani)) +
#   geom_point()
# ggplot(data=dispersionReps, aes(x=rep_num, y=quad, colour=plot_mani)) +
#   geom_point()
# ggplot(data=dispersionReps, aes(x=rep_num, y=yr10, colour=plot_mani)) +
#   geom_point()

ggplot(data=dispersionReps, aes(x=rep_num, y=yr9)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Dispersion Change') +
  scale_x_continuous(breaks=seq(0,50,10)) +
  coord_cartesian(xlim=c(0,45))
#export at 900x900






###look at five factor manipulations for mean change --------------------------------------------------------
#just for the four experiments with five factors, compare to their four factor treatments
# meanFive <- chainsEquations%>%
#   filter(project_name=='e001'|project_name=='e002'|site_code=='NIN'|site_code=='TRA')
#   # filter(treatment=='1_y_n'|treatment=='8_y_n'|treatment=='1_f_u_n'|treatment=='8_f_u_n'|treatment=='2F'|treatment=='3F'|treatment=='4F'|treatment=='ghn'|treatment=='gsn'|treatment=='ncn'|treatment=='nhn'|treatment=='nsn')
# 
# cdr1APlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='A'&variable=='mean'), aes(x=as.factor(treatment), y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
#   # scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#   #                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   # scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(a) CDR e001 A', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr1BPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='B'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   # scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#   #                  labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   # scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(b) CDR e001 B', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr1CPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='C'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   # scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#   #                  labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   # scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(c) CDR e001 C', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr1DPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='A'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   # scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#   #                  labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   # scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(d) CDR e001 D', size=10, hjust='left') +
#   theme(legend.position='none')
# ninPlot <- ggplot(data=subset(meanFive, site_code=='NIN'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('2F', '3F', '4F'),
#                    labels=c('4a', '4b', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(e) NIN herbdiv', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr2APlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='A'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
#   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(f) CDR e002 A', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr2BPlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='B'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(g) CDR e002 B', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr2CPlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='C'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(h) CDR e002 C', size=10, hjust='left') +
#   theme(legend.position='none')
# traPlot <- ggplot(data=subset(meanFive, site_code=='TRA'&variable=='mean'), aes(x=treatment, y=yr9, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('ghn', 'gsn', 'ncn', 'nhn', 'nsn'),
#                    labels=c('4', '4', '4', '5', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'white', 'white', 'black', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(i) TRA lovegrass', size=10, hjust='left') +
#   theme(legend.position='none')
# 
# pushViewport(viewport(layout=grid.layout(2,5)))
# print(cdr1APlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(cdr1BPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(cdr1CPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
# print(cdr1DPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
# print(ninPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 5))
# print(cdr2APlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(cdr2BPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# print(cdr2CPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
# print(traPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
# #export at 2400x1200


#compare any four factor without N to five factor with N
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