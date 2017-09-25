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

rawData <- read.csv('ForBayesianAnalysis_9plusyr_May2017.csv')

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
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp_10\\mv_raw_disp_10_cholesky_3.csv', comment.char='#')
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
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,', '*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=grey) +',
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_9yrplus_greyscale.csv', row.names=F)



###main figure (Figure 1)
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,31,5), labels=seq(1,32,5)) +
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
stat_function(fun=function(x){(-0.4355678805836 + 0.206284709410545*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.646412362975 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.80481178873 + 0.3052417823043*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.881556220806 + 0.28613817272665*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.47239561551705 + 0.169737044522435*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.760383970628 + 0.2131510551534*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.8693507964 + 0.23830501506589*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.64991985234255 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.6200594834885 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.7786483194325 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.9545165927 + 0.184065052587*x + -0.01098335995085*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.79466671152 + 0.221777617765*x + -0.012965195084*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.9937963365 + 0.115927336280465*x + -0.008249695382745*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.749501342645 + 0.1182606338977*x + -0.008947678500279*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-1.0584592386 + 0.0967611136412365*x + -0.0079298397872285*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.82628109123 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.738723922982 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.17290486502852*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.45110273313715 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.787618225855 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.66061839407 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.80453713612 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.721814911125 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6351581713116 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4696952570857 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.60608846617375 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.452870900756935 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6699711315118 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.66157321902375 + 0.20134442028347*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.88592982135 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.53475928527695 + 0.19089452493773*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.418536097163579 + 0.23236697306755*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.91076787185 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.86649543105 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.8970812676 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-1.4090286405 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0.121115667285*x + -0.002017308748187*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0.47146469886428 + 0.177219676375*x + -0.00403474300855*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-0.47322497517325 + 0.3564907277*x + -0.008578629839*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-0.9698554604 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-0.32618030579463 + 0.2767657963*x + -0.006617796859*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-0.62387629564 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0.708044487985 + 0.08193429765005*x + -0.002988068891171*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0.5514306332375 + 0.0928890826627*x + -0.002029429657931*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.002970507965403*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0.37914705462013 + 0.08597626160366*x + -0.0045877395279095*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0.46775451227115 + 0.184512121955*x + -0.0079719394025*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(-0.91609889718 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6190617382075 + 0.18892202129895*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.18554412089625*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.7079711057445 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.47072343405*x + -0.0218738315907*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.60261048651675 + 0.4635459019*x + -0.0227178037995*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3986879669*x + -0.01450388075715*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.46218309836945 + 0.295584306455*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.7757511749545 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6256944194135 + 0.118475996522005*x + -0.00892965167892*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.624271736389 + 0.1277927023669*x + -0.00826777453874*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.746472196975 + 0.1822156758905*x + -0.0106545674046*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.8376894407 + 0.1715947159355*x + -0.0097960821251*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.663151776245 + 0.131526212033795*x + -0.008989751737076*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.5821306307836 + 0.1313266926049*x + -0.008407794926684*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.567791876507 + 0.11309973313418*x + -0.007870093326025*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.8024403945 + 0.179914556634*x + -0.009988968287595*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.73418317324 + 0.1692499436615*x + -0.010184642775198*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.5145035900665 + 0.14405720000105*x + -0.00919920568849*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.7843021853 + 0.10692748972344*x + -0.007068751768034*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.86784299665 + 0.117686228118*x + -0.00810325599218*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.7252108916 + 0.144821980323*x + -0.009614537716765*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.611467512095 + 0.157646301668*x + -0.0099404829764*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-0.80038519715 + 0.10103159627344*x + -0.007330469257587*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(1.51196284075 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(1.61104726795 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(-1.07867602965 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.602502772814982 + 0.32398042698*x + -0.016321782866885*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.31778506881*x + -0.0169378855573065*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.111247854890515*x + -0.00567466532193429*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1368276382228*x + -0.00527377735152355*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.12188872795591*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.102793714747215*x + -0.0046188080437925*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.14461578592426*x + -0.004650863499281*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.134948283307505*x + -0.005324980029944*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.14352495252685*x + -0.007899302774675*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1701261483478*x + -0.007958757626656*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.4896822104435 + 0*x + -0.0055185328601085*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1309599083581*x + -0.00855757322054*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0.09969733873305*x + -0.0054647129962575*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(-1.1410067952 + 0.117466409199155*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(-0.99181219895 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(-0.9711579595 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-1.0392721788 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.44559218644235 + 0.1516238250116*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.857118267975 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.6166930112155 + 0.127066009312165*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.617698133697 + 0.1900112907975*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-1.2593376779 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.0731110144 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.12324908575 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.0827189565 + 0.154443003809913*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.14328092675 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.97039271365 + 0.13143367318762*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.65586066690415 + 0*x + 0.01604005278763*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.42347817691385 + 0*x + 0.012110433940006*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.01664387423458*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.014777369879415*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0132482332508405*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.5650994310945 + 0*x + 0.01638324371943*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.01725144167745*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2412911206505*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.407132602736895 + 0.11961269654436*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.437613436688625 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.153159265995299*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.191787519887195*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.440285739512105 + 0.17021300631072*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.6482900669935 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.5257741723719 + 0.1775519967249*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.50622509000085 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.730932109605 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.62926932949 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.752762155662 + 0.21711896555045*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.9955913878 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.92426765175 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.7054693401903 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.815916808505 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.1500137335 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.750339638542 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.97765712477 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.717407279548 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.22490941585 + 0.10214115133271*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.1644755997 + 0.134896682966898*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.02604375845 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.7993921792 + 0.318668105285*x + -0.01130212019429*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.1839634769 + 0.1598377896485*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-1.19455958745 + 0.1793716280287*x + -0.00660126424872*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.69920793763 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(-1.3642320389 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.38088530265 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.8207061128 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.1082571915 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-0.55985282147765 + 0*x + -0.00554027257574605*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-1.4776208568 + 0.103261432252445*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-0.746240208105 + 0*x + 0*x^2)*(0.1998935)+(0.3711116)}, size=0.5, xlim=c(0,8), colour='grey') +
  #overall lines (average across all datapoints)
  stat_function(fun=function(x){(-0.565723000 + 0.099970750*x + -0.005098725*x^2)*(0.1998935)+(0.3711116)}, size=3, xlim=c(0,31), colour='black')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,31,5), labels=seq(1,32,5)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.4, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.9739324224445 + 0.094844501469285*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-1.8523678852 + 0.1313508039029*x + -0.00374125548844675*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-2.1268778005 + 0.130156815327356*x + -0.00364324991646*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(1.03643245402 + -0.160134718155*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(-0.619708932508599 + 0.189957014099*x + -0.0066855862569*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3039877154*x + -0.0157227177025*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + 0.24900379639*x + -0.0100602111035*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1657606174424*x + -0.00699055407517*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.07541228745 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.35362818165 + -0.233542642591971*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.271260805327583*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.35853447854*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(0 + -0.15099307110809*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + -0.244804681504*x + 0.00785902781389*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.739578604874365 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1760962387693*x + 0.008837530439192*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.767020616542 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.7673977215854 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.85121970351935 + -0.162683048955955*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.71014167888665 + -0.1930306633473*x + 0.00910994574553*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.77970661059275 + -0.16396415511883*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.22799417550505*x + 0.0099497439030025*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.773825033055 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.303544332635*x + 0.0141503962541*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.704122309493295 + -0.203687334919085*x + 0.008793721110144*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.67819730169645 + -0.1896743086813*x + 0.00887513253332085*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.802030633433 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.6534186513918 + -0.23540593282195*x + 0.010442826754635*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.6446535946681 + -0.15106567663135*x + 0.00815000854072045*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.21907102199425*x + 0.01376533283468*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.78979373210835 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.3198425318 + 0.18245715091527*x + -0.018104011062853*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.3287062496 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.1345486883 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.30051826425 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.33989284005 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.3847163156 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.23654579185505*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.184376515690075*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.168924544226315*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2667687933208*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.033084721434105*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.29922661014615*x + 0.0329896630038942*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3505832078526*x + 0.037423322590555*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.9225999830765 + 0*x + 0.05051297978425*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4257919413961*x + 0.059283006314*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3557336172886*x + 0.0559445686252856*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0530064321724415*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0480676288722*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.331256841664585*x + 0.05237257032635*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.182096933056285*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.51784738482485 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.135361958273165*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1000998)+(-0.006758655)}, size=0.5, xlim=c(0,8), colour='grey') +
  #estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)}, size=3, xlim=c(0,31), colour='black')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,31,5), labels=seq(1,32,5)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0.642827022414 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.87928698356595 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.6185675921909 + -0.24223271840166*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(1.1004534008436 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.22713656403045*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.23525048109298*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.21903819673462*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.223809591319153*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.7585067186965 + -0.136546156219445*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.452386649127795 + -0.167146459462145*x + 0.0083193347579525*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.74947611136 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.6819679483915 + -0.12388785588656*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(1.1230206212 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0.84083444856385 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.59463883033195 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.6075560605516 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.051580818355 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.1889691533 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.8795656601925 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.948847717565 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.59981588552745 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.486793758398605 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-0.92765694395 + -0.084129720834135*x + 0.0033636404079775*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-1.12869182705 + -0.09014810296741*x + 0.003614972540545*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + -0.26601444815*x + 0.0066076078085*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(1.22736219455 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + -0.24399921395*x + 0.006019037915*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0.43674323375525 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(-0.82284711479 + -0.1156953942815*x + 0.003620234474879*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(-0.97647088185 + -0.10657681858635*x + 0.0028062735291815*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0.710244378581325 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + -0.08097633172878*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1573503834279*x + 0.004291600883236*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0.66768206427305 + -0.30664170497*x + 0.021463189143885*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.45120533525*x + 0.0265215612222*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.7389647031495 + -0.4737065338*x + 0.02799761667165*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.8351320365355 + -0.2092703095053*x + 0.0161562014652275*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.367313588795*x + 0.021999028384905*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.418702883*x + 0.0239053470931945*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.6266213096652 + -0.34114303503*x + 0.0234959323294*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2897548810865*x + 0.0208890464975121*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.819174005175 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.55182835639325 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.011984049140435*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.42532333344295 + 0.124659075612583*x + -0.0123712904083485*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.58227313193015 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.54184861896185 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.71848630484 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.549570720437 + 0*x + -0.011714261306945*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.444328625148115 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.44186009400385 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.447906112543175 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.42634713764112 + 0.130865311810105*x + -0.0119847037501*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.70546676658 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.5388384280899 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.48609363467948 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.625814115928 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + -0.12141590011795*x + 0.004843743418363*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0.795548808005 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.6943893439175 + -0.269697857266*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.238592906272*x + 0.0121441395529*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.6546792657877 + -0.11226732562184*x + 0.0083046958294845*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.00660862640783635*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.169672529005022*x + 0.010269494615045*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1638955554881*x + 0.0113864860346*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.25903080358*x + 0.017529786834*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2024141212865*x + 0.013153008652*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1808337325653*x + 0.01134805471*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.12087972566032*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17720134319735*x + 0.0095186712054875*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.57260604795 + 0*x + 0.0061543449829295*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.168486481334*x + 0.007240649374781*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0.457944001297835 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.567890671778 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.88950814452 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.85456665005 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.79677369687 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.728906257813 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.21046156508955*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.62498077951 + -0.16083084188332*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.6385767301645 + -0.166652559403575*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2792274268965*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.6440429217785 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2273405658882*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.923413313265 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.649461106545 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.798824228195 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.824933447165 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.924930323975 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.60360001835755 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.023106221818562*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0221860434767331*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0319359244065*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.029453369773905*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03053007068695*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03211800403763*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02576128326638*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.032192290441775*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.45192699929905 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.47584058962112 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.428236151204925 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.6475843255757 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.0189680819 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.595545810256 + 0*x + -0.0083035206967875*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.53692064785725 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.714274762151 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.782476041535 + 0.11150846980336*x + -0.009275295030063*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.1375189225 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.6584192593678 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0.84759935095 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0.63829831215405 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0.4972682333102 + 0*x + 0*x^2)*(0.2505556)+(-0.1346916)}, size=0.5, xlim=c(0,8), colour='grey') +
  #overall line
  stat_function(fun=function(x){(0.401277500 + 0*x + 0*x^2)*0.2505556 + -0.1346916}, size=3, xlim=c(0,31), colour='black')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,31), breaks=seq(0,31,5), labels=seq(1,32,5)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(-0.92037095568735 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.23501812830975*x + 0.01397484562566*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(-0.8145010085785 + -0.27217885264*x + 0.01559658575843*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + -0.18327373893161*x + 0.011846273181816*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + -0.180199259421899*x + 0.01204501985495*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + -0.17868961950286*x + 0.0121915002804055*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6765469699899 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.9759914917505 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7800727730194 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.9170106750475 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0.57064469037927 + 0*x + -0.0030037740657724*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0.984006900105 + 0*x + -0.0035466277264045*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0.194848219905*x + -0.004324598564725*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(-0.70177255018115 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0.20001470909*x + -0.005542491064*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,30), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0 + 0.108694940924144*x + -0.0039088113678813*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0 + 0.11997494370245*x + -0.00431778125233*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,29), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,26), colour='grey') +
  stat_function(fun=function(x){(-0.6841477532121 + 0.300293869856165*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.41141530047*x + -0.02556448191817*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(1.73103614455 + 0.536471298*x + -0.02974354045893*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.769704795437 + 0.20154445349985*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2311908399999*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3495699327*x + -0.015843788873412*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.243672243181645*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.16666570168229*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.15725343954915*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + -0.10896764721473*x + 0.0052702035731655*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(0 + -0.118799791481045*x + 0.005498648533897*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.756305408275 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.9473365727455 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.62122350070505 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1544058798713*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.7801494040755 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.49714530431867 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.835176127136 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.60253151299905 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.59734148247355 + 0.15394085217282*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.63546928433025 + 0.23778137005055*x + -0.013517073470925*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.12442001872101*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.14795755987005*x + -0.009106300140628*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.30729557245*x + -0.01584318014647*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.14372074922026*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1758640570505*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(-0.6186250387779 + 0.140337036374185*x + -0.00852950376224677*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1943440332145*x + -0.0101006911317805*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-0.747074221037 + 0*x + 0*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.118300300146015*x + -0.0076978311035885*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1847339732175*x + -0.00921002796547*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0.15704850426685*x + -0.009270521941685*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(-0.50283811438435 + 0.14051374716395*x + -0.00881607997549785*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0071454555415045*x^2)*(0.1100346)+(0.02923185)}, size=0.5, xlim=c(0,8), colour='grey') +
  #overall mean line
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09597775)+(0.01730774)}, size=3, xlim=c(0,31), colour='black') +
  
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



