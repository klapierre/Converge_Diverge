library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

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

##################################################################################
##################################################################################
#experiment information
expRaw <- read.csv('ExperimentInformation_Mar2016.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('ForBayesianAnalysis_9yr_Nov2016.csv')

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

#raw chains data
chains1 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_0.csv', comment.char='#')
chains1 <- chains1[-1:-334,]
chains2 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_1.csv', comment.char='#')
chains2 <- chains2[-1:-334,]
chains3 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_2.csv', comment.char='#')
chains3 <- chains3[-1:-334,]
chains4 <- read.csv('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\mv_raw_disp\\mv_raw_disp_3.csv', comment.char='#')
chains4 <- chains4[-1:-334,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)

#density plot of chains
plot(density(chainsCommunity$mu_int.1))
plot(density(chainsCommunity$mu_slope.1))
plot(density(chainsCommunity$mu_quad.1))
plot(density(chainsCommunity$U_int.2.1))

#get values for overall (mean) lines across levels of plot mani
#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4)%>%
  gather(key=parameter, value=value, U_int.2.1:mu_quad.4)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))


#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
#set any that are not significant (CI overlaps 0) as 0
chainsIntercept <- chainsCommunity[,8:1159]%>%
  gather(key=parameter, value=value, B.1.1.1:B.4.288.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsIntercept)[1] <- 'parameter'
chainsIntercept <- chainsIntercept%>%
  separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
  mutate(id=as.numeric(id), variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))))%>%
  select(-B, -parameter)

chainsSlope <- chainsCommunity[,1160:2311]%>%
  gather(key=parameter, value=value, B.1.1.2:B.4.288.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsSlope)[1] <- 'parameter'
chainsSlope <- chainsSlope%>%
  separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
  mutate(id=as.numeric(id), variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))))%>%
  select(-B, -parameter)

chainsQuad <- chainsCommunity[,2312:3463]%>%
  gather(key=parameter, value=value, B.1.1.3:B.4.288.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-2*sd, upper=quad+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
  select(parameter, quad)
names(chainsQuad)[1] <- 'parameter'
chainsQuad <- chainsQuad%>%
  separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
  mutate(id=as.numeric(id), variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))))%>%
  select(-B, -parameter)

#merge together with experiment list
chainsExperiment <- chainsIntercept%>%
  left_join(chainsSlope)%>%
  left_join(chainsQuad)%>%
  arrange(id)%>%
  left_join(trtInfo)

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+slope*9+quad*9^2)*(0.1466466)+(0.2925534),
                    ifelse(variable=='dispersion', (intercept+slope*9+quad*9^2)*(0.08732169)+(-0.0001449299),
                           ifelse(variable=='evenness', (intercept+slope*9+quad*9^2)*(0.1008733)+(0.01808835), (intercept+slope*9+quad*9^2)*(0.2157205)+(-0.0550571)))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,', '*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
#write.csv(chainsEquations,'plot mani_equations.csv', row.names=F)



###main figure

#mean change panel
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5808175 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.475868 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4818465 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.276213*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2513305*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.7258905 + 0.332021*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.435645 + 0.337271*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.50645 + 0.346873*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.45086 + 0.3550325*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.333715 + 0.431929*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5741145 + 0.3318175*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.757998 + 0.3328395*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.611494 + 0.328396*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7397655 + 0.4087015*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4880265 + 0.325849*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3380995*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8522285 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.759367 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.10017 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.808065 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6200615 + 0.3250135*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9062515 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.945871 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.683405 + 0.303936*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.85012 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.791925 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.882425 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.73029 + 0.324894*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.858345 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.603426 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4835345 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.345118*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.52635 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.587227 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4890545 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.319503*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.366644*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.622845 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.311657*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3466185*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.651129 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.810761 + 0.406042*x + -0.0543227*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8243225 + 0.420542*x + -0.04928195*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.389315 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.480092 + 0.698731*x + -0.05737565*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.9653115 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.632661*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5829135 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.7912355*x + -0.05771475*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.046595 + 0.806385*x + -0.05473075*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5450195 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.486245*x + -0.03603765*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.918019 + 0.956425*x + -0.0658704*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.543951 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.7447575*x + -0.04288865*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7811385 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5959965 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.819466 + 0.372055*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.12141 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-1.46898 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.334335 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2716195*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3188805*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9267545 + 0.329004*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.369039*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7659355 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.548609 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7325895 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7753005 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.894931 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.21944 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.01564 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.106755 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.154755 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.974013 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5696695 + 0.271251*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4481105 + 0.256148*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4272625 + 0.321772*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4469005 + 0.268282*x + -0.028088*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.485188 + 0.325183*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.521335 + 0.290043*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.368892*x + -0.0316434*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5285185 + 0.247702*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5410335 + 0.235125*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.400328 + 0.3066705*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.649449 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.627494 + 0.211444*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7091685 + 0.2196755*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.58291 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.446766 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.5371295 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(2.139165 + 0.382635*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.80682 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6402225 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.868286 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.40179*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.369857*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.362872*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0.809089 + 0.3469265*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.768032 + 0.409987*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7484105 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.814093 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6353585 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8328295 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.035305 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.129435 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.008565 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9465195 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.874398 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6467465 + 0.415181*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.454267*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.907774 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.025195 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.869424 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00423 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.554558 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2677325*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.511206 + 0.2804*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.472879*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.498339*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3053825*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.340087*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8562845 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7664665 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.748785 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.8234575 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.854959 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.501491 + 0*x + 0.0333599*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.720544 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.460679 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9145405 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.547248 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.08582 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6689675 + 0.3783405*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.743533 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3446785*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.457536*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.8754325 + 0.4724315*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.5067975 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7669195 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4429005*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4287515*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4151295*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.472693*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(-0.885908 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.23181 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7971735 + 0.4439415*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.15661 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.380625 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.106565 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.1987 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.03159 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-1.05086 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.027635 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.05651 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.996121 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9066985 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.488016 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.515935 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6949685 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4787465 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.63018 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6152595 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4668155 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.638643 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7978475 + 0.3720015*x + -0.0517794*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.434865 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.31677 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.942938 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.528255 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1033 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.04882 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.07452 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.079355 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.684575 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.833131 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7231225 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.889898 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.311075*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5441645 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.309692*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9211805 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.531864 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(1.067 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.818075 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7733945 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3438135*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3751525*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6644715 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5871445 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1466466)+(0.2925534)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#last five are the main plot_mani effect lines
  #estimated as mean across treatment lines (plot mani 1-4 staggered by intercept so lines don't overlap)
  #mani1
  stat_function(fun=function(x){(-0.48942950 + 0.20375750*x + -0.01472000*x^2)*(0.1466466)+(0.2925534)}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.47 + 0.20375750*x + -0.01472000*x^2)*(0.1466466)+(0.2925534)}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.45 + 0.20375750*x + -0.01472000*x^2)*(0.1466466)+(0.2925534)}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.43 + 0.20375750*x + -0.01472000*x^2)*(0.1466466)+(0.2925534)}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.48942950 + (0.20375750+0.41946100)*x + -0.01472000*x^2)*(0.1466466)+(0.2925534)}, size=3, xlim=c(0,8), colour='#EC1804')

# print(meanPlot) #export at 1200x1000


#dispersion panel
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.4, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6681845 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.5837085 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.005685 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5288595*x + 0.0706951*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.747229 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8445445 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-2.094725 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.46454*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.649705*x + -0.05429135*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-1.709275 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6661505*x + -0.07759565*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.510603*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.073 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.773341 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.723865 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.276437*x + -0.0445917*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03789665*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0372246*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0343701*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.278634*x + -0.0426509*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-2.42638 + 1.154525*x + -0.110991*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.47081 + 1.06289*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.9873 + 0.994963*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.6994 + 1.07851*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.059 + 1.21855*x + -0.118812*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.469295 + 1.07268*x + -0.1036645*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.61623 + 1.3676*x + -0.137557*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.90899 + 1.264595*x + -0.124802*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.8532645 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8022685 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8620435 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.6782285 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.932494 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.883515 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.643119 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.823011 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.782635 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(-0.8045115 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.30832 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05246945*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.8198095 + 0.722352*x + -0.08747895*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.664444 + 0.5945*x + -0.07640805*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0463844*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4248755*x + -0.0807277*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.886168 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7427205 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7856395 + -0.696456*x + 0.09069225*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.9620225 + -0.5744045*x + 0.06447295*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.673033*x + 0.07835565*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4938675*x + 0.106455*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.619978*x + 0.1109805*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7079895 + -0.748131*x + 0.123813*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.6387025*x + 0.08994005*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.639745 + -0.412311*x + 0.05922375*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.6524695 + -0.4400395*x + 0.0539573*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.7372125 + 0*x + 0.05125795*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.391406*x + 0.0586603*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.364038*x + 0.05478515*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.8632375 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.781623 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.916752 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.904229 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.260195 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7574505 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.803939 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.944868 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.49146 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.44634 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.6599085 + 0.4420585*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.870139 + 0.5010745*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6086905 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.640054 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.13806 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.669601 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.585779*x + 0.0656038*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.519498*x + 0.0824887*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.819161 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.07267 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.6012285*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9114735 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.562854*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.768223 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.02207 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5283155 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.575425 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.28775 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.774136 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.27818 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6581885 + 0*x + 0*x^2)*(0.08732169)+(-0.0001449299)}, size=0.5, xlim=c(0,2), colour='grey') +
#estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)}, size=3, xlim=c(0,8), colour='black')

# print(dispersionPlot) #export at 1200x1000


#richness panel
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6577285 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.647852 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.735546 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.05990855*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.05544875*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.404368*x + 0.06331145*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5787825 + -0.4539835*x + 0.06677535*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.505622 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.542866 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5486295 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5527025 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.408039*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.654859 + -0.2813215*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.017965 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.6153385 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.981623 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.10432 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8227255 + -0.468552*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.63027*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.852865 + -0.430099*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.591159 + -0.399929*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.300345*x + 0.112406*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.9428635 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.884015*x + 0.06779705*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.8970735*x + 0.07751605*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.882181 + -0.598702*x + 0.06650135*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.256185*x + 0.1139675*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.19322*x + 0.1034685*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.951495*x + 0.0898165*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0.8282885 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.702396 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.377171*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.3346765*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.140305 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.6717905 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.639782 + -0.5092635*x + 0.0602285*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0.981696 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(1.02126 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4267305*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.314936*x + -0.0444112*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.408815*x + -0.0491065*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2761505*x + -0.0304935*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2726435*x + -0.03437815*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3139225*x + -0.0389542*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.340755*x + -0.04078155*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.300089*x + -0.04299335*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2231715*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.451024*x + -0.05640925*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.429941*x + -0.0582122*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.444193*x + -0.05933805*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.430044*x + -0.0571693*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3622655*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.42782*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.455888*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7550265 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8490865 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.67758 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.597969 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.304312*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.06173575*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(1.08658 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.119325 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.963229 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.142725 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.724134 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.9991065 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.7520055 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6556145 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.63216 + -0.5119835*x + 0.0664242*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.468817*x + 0.05917585*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.366856*x + 0.0611658*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2880635*x + 0.04205465*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.03898715*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0363922*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2953125*x + 0.0580544*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4821515*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6359855 + -0.4828055*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.6593785 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.8143355 + -0.481316*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.8384195 + -0.4396215*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.36549 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.578715 + 0.5535525*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.7833335 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.59131 + 0.6478105*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.831643 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.7225865*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5956155*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5801155*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8491405 + -0.4202375*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.5588345*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.459266*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4300005*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.986666 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5459325 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.797239 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.578088 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0.6699285 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.691571 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.562259 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.644893 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6646555 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.530787 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4857655*x + 0.07060855*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.5432965*x + 0.0590099*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.892226 + 0.3821975*x + -0.048408*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.798426 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.610528 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.0442279*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.550588 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7256645 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.70977 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.733079 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.722752 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.6775575 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.5663265 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.587569 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.989389 + -0.5353235*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(1.43516 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5180945*x + -0.06500655*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.488131*x + -0.0609247*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3695455*x + -0.05152295*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.336909*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-1.093915 + 0*x + 0*x^2)*(0.2157205)+(-0.0550571)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #plot mani 1-4 staggered to prevent overlap
  #mani1
  stat_function(fun=function(x){(0.23299450 + -0.11865300*x + 0*x^2)*0.2157205 + -0.0550571}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(0.21 + -0.11865300*x + 0*x^2)*0.2157205 + -0.0550571}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(0.19 + -0.11865300*x + 0*x^2)*0.2157205 + -0.0550571}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(0.17 + -0.11865300*x + 0*x^2)*0.2157205 + -0.0550571}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(0.23299450 + (-0.11865300-0.77430550)*x + (0.0092636150+0.06446080)*x^2)*0.2157205 + -0.0550571}, size=3, xlim=c(0,8), colour='#EC1804')

# print(richnessPlot) #export at 1200x1000


#evenness panel
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.578715*x + -0.0559419*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.46274*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0586207*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.053032*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.634606 + -0.3239355*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0425452*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.4965935*x + 0.05976495*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.8290515 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8838255 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.044205 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8414965 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.982405 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.05171805*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.05320535*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0529527*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.05370935*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.435417*x + -0.06163555*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3789635*x + -0.04014615*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.5077505*x + -0.0560322*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.9043835 + 1.677875*x + -0.1353535*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.757219*x + -0.05279425*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.5473125*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7237885 + 0.6809125*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.924575*x + -0.172903*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.944993 + 1.102755*x + -0.08027575*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.8224105*x + -0.0657837*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.657502 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.590741 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8301795 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7925005 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7366995 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.552292 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.34388*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.631626 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.565011 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6169455 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.535477 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.434834*x + -0.04547185*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4279445*x + -0.046662*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.6015905*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4405065*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.6459 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.52507 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.43119 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9815125 + 1.335285*x + -0.09555235*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.76674*x + -0.07029885*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.89171*x + -0.08210835*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.761971*x + -0.0656187*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5970375*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4613355*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.663375*x + -0.0572283*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3267055*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3612005*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.5160325*x + -0.05414345*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.7672245*x + -0.0631065*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.35547*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.474887*x + 0.05600305*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8120335 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.791916 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.761139 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.622859 + 0.361092*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5509265*x + -0.062642*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5990075 + 0.3727385*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6563435 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.71209 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.484256*x + 0.05902065*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4069595*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9864435 + 0.7868795*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9012335 + 0.8165005*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.554648 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1008733)+(0.01808835)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#mean lines by plot mani
  #estimated as mean across treatment lines
  #plot mani 2 and 4 are staggered to prevent overlap
  #mani1
  stat_function(fun=function(x){(-0.21994100 + 0.13346050*x + -0.02002595*x^2)*(0.1008733)+(0.01808835)}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.20 + 0.13346050*x + -0.02002595*x^2)*(0.1008733)+(0.01808835)}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){((-0.21994100+0.32302550) + 0.13346050*x + -0.02002595*x^2)*(0.1008733)+(0.01808835)}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.18 + 0.13346050*x + -0.02002595*x^2)*(0.1008733)+(0.01808835)}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.21994100 + (0.13346050+0.83313400)*x + (-0.02002595-0.05589185)*x^2)*(0.1008733)+(0.01808835)}, size=3, xlim=c(0,8), colour='#EC1804')

# print(evennessPlot) #export at 1200x1000


#print all plots together
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000




###by resource mani
#still need to calculate the proportion of chains where x resource response was greater than y resource response
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

chainsTrt <- chainsEquations%>%
  select(variable, site_code, project_name, community_type, treatment, intercept, slope, quad, yr9, plot_mani, rrich, anpp, MAT, MAP, min_year, experiment_length, alt_length)%>%
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

#plot by resource manipulated at year 9
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='mean'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.20), name='Mean Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(0, 0.6), xlim=c(1,4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.60, label='(a)', size=10, hjust='left')

dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='dispersion'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Dispersion Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(0, 0.15), xlim=c(1,4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.15, label='(b)', size=10, hjust='left')

richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='richness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.25, 0.2), xlim=c(1,4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.2, label='(c)', size=10, hjust='left')

evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='evenness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr9', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Evenness Change') +
  scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
                   labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.06, 0.06), xlim=c(1,4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.06, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600


      
###summary stats from bayesian output
#gather summary stats needed and relabel them
chainsCommunitySummary <- chainsCommunity%>%
  select(U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4, E_int.1.1, E_int.2.1, E_int.3.1, E_int.4.1, E_slope.1.1, E_slope.2.1, E_slope.3.1, E_slope.4.1, E_quad.1.1, E_quad.2.1, E_quad.3.1, E_quad.4.1, E_int.1.2, E_int.2.2, E_int.3.2, E_int.4.2, E_slope.1.2, E_slope.2.2, E_slope.3.2, E_slope.4.2, E_quad.1.2, E_quad.2.2, E_quad.3.2, E_quad.4.2, D_int.1.1, D_int.2.1, D_int.3.1, D_int.4.1, D_slope.1.1, D_slope.2.1, D_slope.3.1, D_slope.4.1, D_quad.1.1, D_quad.2.1, D_quad.3.1, D_quad.4.1, D_int.1.2, D_int.2.2, D_int.3.2, D_int.4.2, D_slope.1.2, D_slope.2.2, D_slope.3.2, D_slope.4.2, D_quad.1.2, D_quad.2.2, D_quad.3.2, D_quad.4.2)%>%
  gather(key=key, value=value, U_int.2.1:D_quad.4.2)%>%
  group_by(key)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(CI=sd*2)%>%
  separate(key, into=c('head', 'variable', 'level'), sep='\\.', remove=F)%>%
  separate(head, into=c('type', 'parameter'))%>%
  mutate(variable=ifelse(variable==1, 'mean change', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))))%>%
  mutate(type_level=paste(type, level, sep='_'))%>%
  mutate(predictor=ifelse(type_level=='E_1', 'ANPP', ifelse(type_level=='E_2', 'gamma diversity', ifelse(type_level=='D_1', 'MAP', ifelse(type_level=='D_2', 'MAT', ifelse(type_level=='U_1', 'plot mani 2', ifelse(type_level=='U_2', 'plot mani 3', ifelse(type_level=='U_3', 'plot mani 4', ifelse(type_level=='U_4','plot mani 5', 'overall')))))))))%>%
  select(parameter, variable, predictor, median, sd, CI)

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)


#mean plots
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='int'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='slope'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quad'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.1, 0.1) +
  coord_flip()

#dispersion plots
dispersionIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='int'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

dispersionSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='slope'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

dispersionQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quad'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.1, 0.1) +
  coord_flip()

#richness plots
richnessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='int'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

richnessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='slope'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

richnessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quad'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.1, 0.1) +
  coord_flip()

#evenness plots
evennessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='int'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

evennessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='slope'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

evennessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quad'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma diversity', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.1, 0.1) +
  coord_flip()

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


#overall responses
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean change' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('\nMean Change') +
  annotate('text', x=3.45, y=-0.8, label='(a)', size=10, hjust='left')

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='dispersion' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.15), breaks=seq(-0.2, 0.2, 0.2)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
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
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
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
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.4, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
#export at 2400x500





###look for patterns of spp appearance/disappearance -- no clear patterns, probably because just the few CDR examples that are long term enough to see the pattern
relAbund <- read.csv('SpeciesRelativeAbundance_Nov2016.csv')%>%
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

ggplot(data=turnoverDiff, aes(x=as.factor(plot_mani), y=turnover)) +
  geom_boxplot() +
  xlab('Number of Factors Manipulated') +
  ylab('Species Turnover')
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





#look at number replicates for dispersion results -- doesn't make a difference
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
  




###look at five factor manipulations for mean change
# #just for the four experiments with five factors, compare to their four factor treatments
# meanFive <- mean4%>%
#   filter(treatment=='1_y_n'|treatment=='8_y_n'|treatment=='1_f_u_n'|treatment=='8_f_u_n'|treatment=='2F'|treatment=='3F'|treatment=='4F'|treatment=='ghn'|treatment=='gsn'|treatment=='ncn'|treatment=='nhn'|treatment=='nsn')
# 
# cdr1APlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='A'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
#   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#                      labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(a) CDR e001 A', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr1BPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='B'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(b) CDR e001 A', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr1CPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='C'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(c) CDR e001 A', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr1DPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='A'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(d) CDR e001 D', size=10, hjust='left') +
#   theme(legend.position='none')
# ninPlot <- ggplot(data=subset(meanFive, site_code=='NIN'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('2F', '3F', '4F'),
#                    labels=c('4a', '4b', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(e) NIN herbdiv', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr2APlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='A'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
#   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(f) CDR e002 A', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr2BPlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='B'), aes(x=treatment, y=yr10, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(g) CDR e002 B', size=10, hjust='left') +
#   theme(legend.position='none')
# cdr2CPlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='C'), aes(x=treatment, y=final_year_estimate, fill=treatment)) +
#   geom_bar(stat="identity", colour='black') +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
#                    labels=c('4', '5')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'black')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(h) CDR e002 C', size=10, hjust='left') +
#   theme(legend.position='none')
# traPlot <- ggplot(data=subset(meanFive, site_code=='TRA'), aes(x=treatment, y=final_year_estimate, fill=treatment)) +
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

meanCompare <- subset(resourceMani, variable=='mean')%>%
  select(variable, site_code, project_name, community_type, treatment, plot_mani, yr9)%>%
  left_join(expRawMean, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))


#plot without N at four factors, with N at five factors
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none')
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&herb_removal>0), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none')
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plant_mani>0), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
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