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

rawData <- read.csv('ForBayesianAnalysis_9-20 year_May2017.csv')

rawData2<- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  summarise(mean_mean=mean(mean_change), std_mean=sd(mean_change), mean_disp=mean(dispersion_change), std_disp=sd(dispersion_change), mean_rich=mean(S_PC), std_rich=sd(S_PC), mean_even=mean(SimpEven_change), std_even=sd(SimpEven_change)) #to backtransform

#select just data in this analysis
expInfo2 <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani))

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
# chains1 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_9-20\\mv_raw_disp_9-20_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_9-20\\mv_raw_disp_9-20_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_9-20\\mv_raw_disp_9-20_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_9-20\\mv_raw_disp_9-20_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)
# 
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsCommunity$mu.1.1))
# plot(density(chainsCommunity$mu.1.2))
# plot(density(chainsCommunity$mu.1.3))


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
# write.csv(chainsCommunity2, 'bayesian_output_summary_9yrplus_09252017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_9yrplus_09252017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,4648:6495]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 4648:6495])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,4648:6495]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 4648:6495])'] <- 'sd'
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
# write.csv(chainsFinal, 'bayesian_output_mean sd_9yrplus_09252017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_9yrplus_09252017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=20, 19, alt_length))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,', '*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=grey) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_9yrplus_greyscale.csv', row.names=F)



###main figure (Figure 1)
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(0,19,5), labels=seq(1,20,5)) +
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(-0.435581514294425 + 0.21916367106553*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.587363872388 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.778000517914 + 0.3181098058085*x + -0.0146451034951525*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.92824179521 + 0.314945226716*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.4984123146268 + 0.20905793944106*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.7885912412045 + 0.2720028073842*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.88548114437 + 0.309109082062*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.774857333765 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.737562332668 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.918177457305 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.96697435555 + 0.19793772213965*x + -0.011704651938895*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.84251092255 + 0.263510926355*x + -0.0158964427121*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.9846962015 + 0.122818829097195*x + -0.00875300259366265*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.741506648108 + 0.137464210627808*x + -0.010528042624265*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-1.0017139283 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.8091782716085 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.696392430532375 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.751639407515 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.66531570881265 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.787063657315 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.695901647295 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.58518650251655 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.42216733906315 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.554985313508025 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.43115274079095 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6396019477765 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.647201269665 + 0.21542716919745*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.819405301165 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.491586765902885 + 0.20523783401595*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.26828298676225*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.877337830905 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.827862499775 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.8767744374 + 0.122935531090005*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-1.37506769 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.23889903707*x + -0.0081320614920694*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.265321483785*x + -0.00795795516121985*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.8830462036 + 0.53190397405*x + -0.017497455071*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-1.17492461495 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.96597272255 + 0.5449302663*x + -0.020932278055*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.5996339997325 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.46652532796313 + 0.218228085445*x + -0.0106740996707*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.237169668595*x + -0.0105604644415*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.105446831284105*x + -0.0054724084435345*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.28115147735*x + -0.0162691844415*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.30043388735*x + -0.013892972002*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.856185576435 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.7300407095495 + 0.25706121940595*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.415916944029895 + 0.239200020101*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.598713545974454 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.60033213275*x + -0.0326635143145*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.46447256139695 + 0.606466443*x + -0.035529399354*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.49506391385*x + -0.022854451705775*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6133762060306 + 0.38563704205*x + -0.017984295011255*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.81376274983 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.58707481851045 + 0.111198713326165*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.6047827555605 + 0.13629152982881*x + -0.008906845929477*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.749149523745 + 0.211787067382*x + -0.013103953292265*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.81168044825 + 0.1780538228919*x + -0.010248897862598*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.599433227839 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.57533288328085 + 0.144566406007231*x + -0.009404306962934*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.55239052789522 + 0.120309262105085*x + -0.008292289171944*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.80008663879 + 0.19185449255285*x + -0.010592043581164*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.721368755065 + 0.18011221677725*x + -0.010905639204365*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.5108607931855 + 0.172544923505605*x + -0.01165366143294*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.7610259979 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.79621348145 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.7108810440725 + 0.14995991418565*x + -0.009808884812018*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.59385327085 + 0.1654352921691*x + -0.010317786144494*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.78003153935 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(1.5708956372 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(1.59260583285 + 0.141241132113095*x + -0.006354647969635*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.9838440487 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.611040570717 + 0.357211965095*x + -0.018943532993696*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.355330706565*x + -0.0199075939272531*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.18907399233885*x + -0.01093751027566*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.26947242088*x + -0.0153514009934*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.19717092789049*x + -0.010143758624179*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.190593421988335*x + -0.011238496367057*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.26945644518*x + -0.01401743428095*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2930837927595*x + -0.0168648122712*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.295733029589*x + -0.0182671397775*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.33731728865*x + -0.019370298882*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1641278858998*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.15475614435163*x + -0.01187150137293*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2117541068763*x + -0.015097169957445*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31536469004*x + -0.0219401670465*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.20547024327075*x + -0.013097668730587*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-1.16248991165 + 0.130182803671*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(-0.956910693 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(-0.96566803705 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-1.03211731485 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.48462956621955 + 0.20490065653785*x + -0.011967075090137*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0.14237285533005*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.84175343291 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.7107839771935 + 0.1892863849217*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.6839077142475 + 0.261088089175*x + -0.0129600397416785*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-1.232194437 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.05332812625 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.13752295055 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.1345952019 + 0.18402222253425*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.13998517955 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.99087283265 + 0.15615761192452*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.61141129608735 + 0*x + 0.017962239438545*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.4429812585293 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0180176201159625*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.016553998814273*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.4950106805827 + -0.141614655409006*x + 0.018499683950605*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0201421271859*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.274873194789*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.38340022630215 + 0.12840886150752*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.397600031083485 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.175402249309755*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2365431736735*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.436450277166845 + 0.19121854186099*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5530496034625 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.526481631087 + 0.198395785376792*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.44691747968345 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7142815184975 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5792829111356 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7672893045566 + 0.241623584942*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.9337879778 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.9108426283 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6791805191225 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.77428246312535 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.1030920047 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6992874054945 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.03645024345 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6283362339695 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.2277181416 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.1700126941 + 0.145405212676535*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.00753995595 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.931877886 + 0.39875464119*x + -0.01669923436634*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.22008148585 + 0.1881240265944*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.2347192331 + 0.20736551983145*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.6422861990415 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.33969377485 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.35261825015 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.8029787937 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.98656205925 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-0.47771079508751 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-1.6295964554 + 0.14500952494331*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-0.6664927361045 + 0*x + 0*x^2)*(0.1930012)+(0.3613757)}, size=0.5, xlim=c(0,8), colour='grey') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.560286000 + 0.112800500*x + -0.005951565*x^2)*(0.1930012)+(0.3613757)}, size=3, xlim=c(0,19), colour='black')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(0,19,5), labels=seq(1,20,5)) +
  # scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.4, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.72364616190335 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-1.26079542925 + -0.16115458767275*x + 0.01379813838125*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-1.56430904685 + -0.20932186114705*x + 0.0180558599414*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(1.0295207950175 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.159242201151175*x + 0.009209219377486*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.26359256572*x + 0.013277408321235*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.30907849286*x + -0.0153286294199*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.21455782431395*x + -0.0123049427171915*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.30272109719*x + -0.01314716586225*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.27217830908*x + -0.0144392403114*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.8594142423805 + -0.31191467399265*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.1229912587 + -0.402545897298*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.260369035146615*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.429744432511*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.38063373321125*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1637102793844*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.70396704274529 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2424516321065*x + 0.013295103141615*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.7913718719525 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.77809469236085 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.79164035725255 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.72192897295208 + -0.227464808493385*x + 0.0113185095668215*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.752576709595 + -0.208941991794265*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.53614857107345 + -0.33467992085015*x + 0.01831286461932*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.7590837082175 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.60026118709921 + -0.38134710082*x + 0.01908298280665*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.76453268031765 + -0.271098889357*x + 0.0132601139918415*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.7577815594145 + -0.25429895963855*x + 0.013068806294715*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.84539709076805 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.63595246276955 + -0.2763267003799*x + 0.014049037782861*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.64862975297265 + -0.20263801263016*x + 0.0125187974938665*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2578665812375*x + 0.015962802753425*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.776739937943 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.23466108415 + 0.24909237928206*x + -0.024275982148515*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.2664332109 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.06453565309 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.2504230563 + 0*x + -0.0178625486313065*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.24799598465 + 0*x + -0.0182676037955385*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.32360156155 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.25397995427795*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.202799663514785*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.304953319448*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0393810397225375*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.04550668549955*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.347276427650455*x + 0.0457878220858*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.672027341708885 + -0.40529162471265*x + 0.046381227610385*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.470385582135*x + 0.0518491168661*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.97384749162725 + -0.323716540703955*x + 0.058915175578105*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.536245884578*x + 0.072136976775*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4456990220535*x + 0.06728509484*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3937742230877*x + 0.066016292575*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.05564645958413*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.419549182015*x + 0.0627827693116*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.382035715782395*x + 0.049858415854895*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.20113390838195*x + 0.01977969401641*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2425233617177*x + 0.022387687616687*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2588813355619*x + 0.0225176871875065*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.62271354134565 + 0*x + 0*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.49268403737225 + -0.1921556639636*x + 0.0186082291106*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3194301961945*x + 0.02695264559583*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0185028213156535*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1906747945454*x + 0.020915364073275*x^2)*(0.09793173)+(-0.004155456)}, size=0.5, xlim=c(0,8), colour='grey') +
  #estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + -0.107810500*x + 0.014672050*x^2)*0.09793173-0.004155456}, size=3, xlim=c(0,19), colour='black')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(0,19,5), labels=seq(1,20,5)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.72862618389195 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.8319769269585 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.59334127812505 + -0.250297795091285*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.979969583265 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.243447361898908*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.319243840896655*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.32486654384985*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.24703985818447*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.837578655555 + -0.17008138868105*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.54485632867505 + -0.2157152845375*x + 0.0116200293696253*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.7265167632765 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.670371536087 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(1.12548237855 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.8164978983695 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.569237451894495 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.49225599433759 + 0.26257136641035*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.11676106805 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.23522020235 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.956252795755 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.512743863656445 + -0.305661341190062*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.00848228775 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.605650238441815 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.47329729462395 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.481830479505835 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.011484410460763*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0.7117107339265 + -0.31844820715*x + 0.0171210147915*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.57094713805*x + 0.031317678705*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5792103084*x + 0.032841406295*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.392385694*x + 0.0136211345084*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(1.43099590985 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.98863332535 + -0.4447963989*x + 0.016060631332*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.674467791652 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4072008748*x + 0.0199606210675*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4112092767*x + 0.019949001408*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.4754523718146 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17084135316625*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.6173948519025 + -0.3329881921105*x + 0.0280741647104*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6730155814*x + 0.04568089922*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.7143656626*x + 0.0488439813765*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.7308372340305 + -0.210729997329005*x + 0.02082265711336*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5620081482*x + 0.0391244989775*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.66862948995*x + 0.04425683973*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.561459747933855 + -0.4710185353*x + 0.036248719138*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4590881913*x + 0.03605855937*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.223615489161855*x + 0.024867091144535*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.79746761984 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.51159435246495 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1390957255242*x + -0.015262344551819*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0.16949348666987*x + -0.016006531530315*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.587921382 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.52197538891605 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.723797984046 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.529892999012015 + 0*x + -0.01447536543217*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.398670926834 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.4335771149224 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0.17333425634374*x + -0.0153953959390745*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.70337984239 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.563294116287445 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.4542461215626 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.733386280337 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.4493182436457 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2080215007645*x + 0.01023206048223*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.70254110428 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.7534084139255 + -0.2995049293399*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.283087258401*x + 0.014055157910305*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.56289391020557 + 0*x + 0.0116833309831055*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.24183289324035*x + 0.01573188778287*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.161751822448926*x + 0.0110491881395655*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.33477263489*x + 0.0206801441155*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3398626104935*x + 0.02329065885*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.180426445763865*x + 0.0103757889880597*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17356207102557*x + 0.0098678784310245*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0.4908685780489 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.19424975569245*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0.434221440952425 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.5229171504573 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.46563718554685 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.944899625075 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.92875838915 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7751114318945 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.72754131481 + -0.169284673356655*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.27820160837215*x + 0.017216628078955*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.5776364439455 + -0.16695713226981*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.6630673178598 + -0.2131345472535*x + 0.0144146738999935*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.354642815605*x + 0.019026542754569*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.6128343193688 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2788241038325*x + 0.018272855197935*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.7786298281505 + 0*x + -0.01685205056015*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.45594889457755 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.50880386100355 + 0*x + -0.0146222905642851*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.7296076320215 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.774310450298 + 0*x + -0.0141840775815285*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.8213167183525 + 0*x + -0.0159483092703674*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.542094750186645 + 0*x + 0.0227186575459235*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0237303469400393*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0240973165162074*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.026219854830202*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3081700397896*x + 0.0463614672732*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.241070952034945*x + 0.037661686642165*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.26806314766*x + 0.0397808543553*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.22957629180127*x + 0.0394076604391*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03358836696819*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0319026853261935*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0331523487099705*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.80044430955755 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.4099668743062 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0138678232146071*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.58056776715705 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.0262684844 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.49011995543119 + 0.16974297114314*x + -0.0150651922567005*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.45788435791919 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.61192858510375 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.70317032842985 + 0.18879802866815*x + -0.0174211064788675*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.07549810525 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.6044416304095 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0.753542377085 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0.5061099608946 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0.4454577144888 + 0*x + 0*x^2)*(0.243293)+(-0.1275627)}, size=0.5, xlim=c(0,8), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.436278500 + -0.079491550*x + 0*x^2)*0.243293 + -0.1275627}, size=3, xlim=c(0,19), colour='black')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,19), breaks=seq(0,19,5), labels=seq(1,20,5)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.76456427484555 + -0.31467876011736*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.256298813377*x + 0.01562113446646*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.6316573532558 + -0.408463245115*x + 0.0261247372895*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1870834522009*x + 0.012070056291387*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.432910027272685 + -0.20656105667665*x + 0.014437552587645*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.567263864730945 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.656954970857084 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.8283367294355 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7145065495872 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.8108281774765 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2428179476487*x + -0.018702920434182*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0.24702617177885*x + -0.02021363362557*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.4866251759268 + 0.2259064889476*x + -0.0104975569954405*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.44790311165*x + -0.025532134025*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.65910802635*x + -0.04141798077*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.29463035334*x + -0.010358604735778*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.730860213387 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.20509497837675*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.475215196637055 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(-0.46694053112985 + 0.263190405295*x + -0.01107841448891*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.290281866765*x + -0.0123578044653065*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1360289837169*x + -0.0071922644592477*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.65208497645*x + -0.04868567786525*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.82047266675425 + 1.0888206242*x + -0.07927745034*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6490281450695 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6096691747449 + 0.3444174915266*x + -0.023514006658357*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.64347494925*x + -0.039563693028*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3585799194605*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.5487485781418 + 0.2145957582389*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.5532888451756 + 0.19322894575057*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + -0.18624188862218*x + 0.010369002031167*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.22991272392685*x + 0.01343114727695*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1921001108446*x + 0.01085321066214*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6730527026149 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.857531422687 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.230134760619531*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17853132700551*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,19), colour='grey') +
  stat_function(fun=function(x){(0 + 0.149218352823705*x + -0.00851289417905*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0.21023651193913*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2254157330236*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.210947975825515*x + -0.019364900029615*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.49474211330083 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.6092940485594 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.43913926471636 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.6400389870055 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.4907774582232 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.492528801724885 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5226810959842 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.52845506943009 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5503199950775 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.44255278294*x + -0.026112151202117*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.219395214109895*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.4433177112303 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.19094669532795*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.5474998992231 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0.26608275942436*x + -0.0176613916831315*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=0.5, xlim=c(0,8), colour='grey') +
  #overall mean line
  stat_function(fun=function(x){(-0.183515500 + 0*x + 0*x^2)*(0.1087544)+(0.02925361)}, size=3, xlim=c(0,19), colour='black')
  
  # print(evennessPlot) #export at 1200x1000
  
  
  #print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000


# ##summary stats from bayesian output --------------------------------------------------------
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
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_9yrplus_09252017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_9yrplus_09252017.csv')

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
  scale_y_continuous(limits=c(-1, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=40, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Mean\nChange') +
  annotate('text', x=3.45, y=-1, label='(a)', size=10, hjust='left')

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='dispersion' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.4), breaks=seq(-0.2, 0.2, 0.2)) +
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
  scale_y_continuous(limits=c(-0.3, 0.7), breaks=seq(-0.3, 0.5, 0.3)) +
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



