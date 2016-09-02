library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

setwd('C:\\Users\\Kim\\Desktop\\bayesian output')

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
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), precip=mean(precip))

rawData <- read.csv('9 yr subset\\ForBayesianAnalysis_9yr_abs value_Aug2016.csv')

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
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich), anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  summarise(length_median=median(experiment_length), length_min=min(experiment_length), length_max=max(experiment_length),
            plot_mani_median=median(plot_mani), plot_mani_min=min(plot_mani), plot_mani_max=max(plot_mani),
            rrich_median=median(rrich), rrich_min=min(rrich), rrich_max=max(rrich),
            anpp_median=median(anpp), anpp_min=min(anpp), anpp_max=max(anpp),
            MAP_median=median(MAP), MAP_min=min(MAP), MAP_max=max(MAP),
            MAT_median=median(MAT), MAT_min=min(MAT), MAT_max=max(MAT)
            )%>%
  gather(variable, estimate)

################################################################################
################################################################################

#raw chains data
chains1 <- read.csv('9 yr subset\\diversity_absvalue_MV_rdisp_9yr_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('9 yr subset\\diversity_absvalue_MV_rdisp_9yr_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('9 yr subset\\diversity_absvalue_MV_rdisp_9yr_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('9 yr subset\\diversity_absvalue_MV_rdisp_9yr_3.csv', comment.char='#')
chains4 <- chains4[-1:-5000,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)

#get values for overall (mean) lines across levels of plot mani
#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4)%>%
  gather(key=parameter, value=value, U_int.2.1:mu_quad.4)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))



###mean change
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsMeanIntercept <- chainsCommunity[,7:295]%>%
  gather(key=parameter, value=value, B.1.1.1:B.289.1.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsMeanIntercept)[1] <- 'parameter1'

chainsMeanSlope <- chainsCommunity[,1163:1462]%>%
  gather(key=parameter, value=value, B.1.1.2:B.289.1.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsMeanSlope)[1] <- 'parameter2'

chainsMeanQuad <- chainsCommunity[,2319:2608]%>%
  gather(key=parameter, value=value, B.1.1.3:B.289.1.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-2*sd, upper=quad+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
  select(parameter, quad)
names(chainsMeanQuad)[1] <- 'parameter3'


#merge together with experiment list
mean <- cbind(chainsMeanIntercept, chainsMeanSlope, chainsMeanQuad)%>%
  select(-parameter2, -parameter3)%>%
  separate(parameter1, into=c('head', 'exp', 'type', 'variable'), sep='\\.', remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

mean2 <- mean[order(mean$exp2),]
mean3 <- cbind(expInfo2, mean2)%>%
  select(-exp2, -exp)
  
mean4 <- left_join(mean3, expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  #get estimates at various time points
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.1462826+0.291631,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.1462826+0.291631,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.1462826+0.291631)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(mean4,'mean_equations.csv', row.names=F)



#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.634485 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3012365*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3107295*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.787923 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.41016 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.5389 + 0.2974235*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.41791 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.4253 + 0.341261*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5903105 + 0.3155305*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.828467 + 0.2788105*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.6803725 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3734695*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.724843 + 0.3851065*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.486189 + 0.3147705*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.309967*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.870335 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.726106 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.0961 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8191465 + 0.237924*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.62252 + 0.3229085*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.871017 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.879785 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.664205 + 0.2845315*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.814075 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.716955 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.7896 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.717535 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.81288 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.618591 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2888875*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.563536 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.530707 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5482295 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.310588*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2741825*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5974145 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.276746*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3377865*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.613955 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7676155 + 0.325834*x + -0.0422435*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7756375 + 0.352153*x + -0.03953805*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.29719 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0.5896945 + 0.572547*x + -0.04264395*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-1.0135 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.5805545 + 0.7338565*x + -0.0325476*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.593862 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.7531585*x + -0.05217525*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.169625 + 0.7416355*x + -0.04285035*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4316245*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.10302 + 0.8560295*x + -0.05613695*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5639765 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.738684*x + -0.04006805*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7746965 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5789545 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.717084 + 0.3323745*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.090025 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-1.37877 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.255225 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3719675*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.897104 + 0.363031*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3889695*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6922085 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5997145 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.746746 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.761638 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.868855 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.131805 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9707205 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.075995 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.12862 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9724805 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5361945 + 0.2505775*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2991405*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.238461*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.52457 + 0.330101*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.474252 + 0.278784*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.369088*x + -0.0293092*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5060315 + 0.237116*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4395195 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2846355*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5719695 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.521159 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.669847 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.281308*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2832025*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6534675 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(2.20815 + 0.311397*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.808552 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7379305 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8789 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.434565*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3476265*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.342158*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0.924302 + 0.282543*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.784 + 0.3938115*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.747369 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8244855 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.630483 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.748254 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00324 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.08183 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9507325 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8708825 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.774347 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.642302 + 0.401849*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4302045*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.906543 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.09382 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8596645 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.046065 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.558444 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2566905*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.5537035 + 0.2683935*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4436315*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4471005*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3096775*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.344421*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8589365 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8011195 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7514715 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.785334 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.542423 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.464873 + 0.2440575*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.859473 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.576093 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.839049 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.496052 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.989648 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.57072 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9489815 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.549743 + 0.2946155*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.723641 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.751209 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3453985*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5414045 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.504778 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.886894 + 0.328054*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.388722*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.508559 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.758133 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.405138*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3529605*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3611935*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3978935*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7811395 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.10908 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.657341 + 0.388154*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.128405 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.34123 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.11346 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.180905 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.0364 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-0.974548 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.94649 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9718855 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9807975 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8908295 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.500475 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.647234 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6334475 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5999145 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4759365 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.629921 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8269555 + 0.3103255*x + -0.04019265*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.479655 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.368745 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.921859 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.00872271*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.71178 + 0.252229*x + -0.0141383*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.005825 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.04072 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9983215 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.04534 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.596306 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7452885 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.746615 + 0.311881*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.035195 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.538404 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2993405*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.883181 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4865685 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.985571 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.755507 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.699147 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.371415*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4238005*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.95118 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.610022 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0.6862835 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1462826)+(0.291631)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
#last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.4833985 + 0.1832865*x + -0.0119574*x^2)*0.1462826 + 0.291631}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.4833985 + 0.1832865*x + -0.0119574*x^2)*0.1462826 + 0.291631}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.4833985 + 0.1832865*x + -0.0119574*x^2)*0.1462826 + 0.291631}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.4833985 + 0.1832865*x + -0.0119574*x^2)*0.1462826 + 0.291631}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.4833985+0.6136355) + (0.1832865+0.4156465)*x + -0.0119574*x^2)*0.1462826 + 0.291631}, size=3, xlim=c(0,8), colour='#EC1804')

# print(meanPlot) #export at 1200x1000




###dispersion
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsDispersionIntercept <- chainsCommunity[,296:584]%>%
  gather(key=parameter, value=value, B.1.2.1:B.289.2.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsDispersionIntercept)[1] <- 'parameter1'

chainsDispersionSlope <- chainsCommunity[,1452:1750]%>%
  gather(key=parameter, value=value, B.1.2.2:B.289.2.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsDispersionSlope)[1] <- 'parameter2'

chainsDispersionQuad <- chainsCommunity[,2608:2896]%>%
  gather(key=parameter, value=value, B.1.2.3:B.289.2.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-2*sd, upper=quad+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
  select(parameter, quad)
names(chainsDispersionQuad)[1] <- 'parameter3'


#merge together with experiment list
dispersion <- cbind(chainsDispersionIntercept, chainsDispersionSlope, chainsDispersionQuad)%>%
  select(-parameter2, -parameter3)%>%
  separate(parameter1, into=c('head', 'exp', 'type', 'variable'), sep='\\.', remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

dispersion2 <- dispersion[order(dispersion$exp2),]
dispersion3 <- cbind(expInfo2, dispersion2)%>%
  select(-exp2, -exp)

dispersion4 <- left_join(dispersion3, expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  #get estimates at various time points
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.06185618+0.06370134,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.06185618+0.06370134,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.06185618+0.06370134)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(dispersion4,'dispersion_equations.csv', row.names=F)



#main figure
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(0,0.5))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.5, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7098505 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.727042 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.675392 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.680259 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8157985 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7380995 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7723565 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.751919 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.658329 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.39195 + 0.5676825*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3890605*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.282115 + 0.5073785*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.532906*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-0.699862 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0.9631295 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(2.179385 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.550475 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.21146 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.8263025 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(1.221125 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.557935 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.9491775 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.703444 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.6126005*x + -0.0686756*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0.772649 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.747577 + 0.368461*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7862855 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.047949*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.737552 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.778986 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.3768 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.63489 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4867975*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4121975*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4115385*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.670543 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.698251 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0125573*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.440911*x + -0.0249044*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8398535 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.0181 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.9663715 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(2.35909 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.864375 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.06185618)+(0.06370134)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
#estimated as mean across treatment lines
  #overall lines 
  stat_function(fun=function(x){(-0.3500895 + 0.1081475*x + 0*x^2)*0.06185618 + 0.06370134}, size=3, xlim=c(0,8), colour='#1400E5') +
  stat_function(fun=function(x){((-0.3500895+0.2349185) + 0.1081475*x + 0*x^2)*0.06185618 + 0.06370134}, size=3, xlim=c(0,8), colour='#4A06AC') +
  stat_function(fun=function(x){(-0.3500895 + 0.1081475*x + 0*x^2)*0.06185618 + 0.06370134}, size=3, xlim=c(0,8), colour='#800C74') +
  stat_function(fun=function(x){(-0.3500895 + 0.1081475*x + 0*x^2)*0.06185618 + 0.06370134}, size=3, xlim=c(0,8), colour='#B6123C') +
  stat_function(fun=function(x){((-0.3500895+0.8597060) + 0.1081475*x + 0*x^2)*0.06185618 + 0.06370134}, size=3, xlim=c(0,8), colour='#EC1804')

# print(dispersionPlot) #export at 1200x1000







###richness
chainsRichnessIntercept <- chainsCommunity[,874:1163]%>%
  gather(key=parameter, value=value, B.1.4.1:B.289.4.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsRichnessIntercept)[1] <- 'parameter1'

chainsRichnessSlope <- chainsCommunity[,2030:2319]%>%
  gather(key=parameter, value=value, B.1.4.2:B.289.4.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsRichnessSlope)[1] <- 'parameter2'

chainsRichnessQuad <- chainsCommunity[,3186:3474]%>%
  gather(key=parameter, value=value, B.1.4.3:B.289.4.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-2*sd, upper=quad+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
  select(parameter, quad)
names(chainsRichnessQuad)[1] <- 'parameter3'


#merge together with experiment list
richness <- cbind(chainsRichnessIntercept, chainsRichnessSlope, chainsRichnessQuad)%>%
  select(-parameter2, -parameter3)%>%
  separate(parameter1, into=c('head', 'exp', 'type', 'variable'), sep='\\.', remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

richness2 <- richness[order(richness$exp2),]
richness3 <- cbind(expInfo2, richness2)%>%
  select(-exp2, -exp)

richness4 <- left_join(richness3, expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  #get estimates at various time points
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.1555338+0.1586638,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.1555338+0.1586638,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.1555338+0.1586638)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(richness4,'richness_equations.csv', row.names=F)



#main figure
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(0,0.8))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.576157 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.595845 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.555378 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5653345 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5799415 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.436655*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5980845 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4727495*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4734655*x + -0.05380915*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4984125*x + -0.05702875*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.569302 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.604432 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.56818 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6067385 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.582312 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.523927 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.542593 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4546985*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.373304*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5668135 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.543287 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.335245*x + -0.105977*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.8605215*x + -0.07257455*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.058035*x + -0.0897439*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.800237*x + -0.0838634*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.364115*x + -0.1157485*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.187205*x + -0.09676145*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.08107*x + -0.09648895*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.430776*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4214415*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.465644*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4699225*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5316335 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.334416*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.630722 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5727405 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.388774*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4277925*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3003385*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3212035*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4078465*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.51564*x + -0.0665281*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4325875*x + -0.05457*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.486367*x + -0.0552819*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.94406 + 0.557373*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.757052 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.873719 + 0.752108*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.8964475*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3070815*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2951935*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5606365 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3922145*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3462205*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4270185*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0.7375435 + 0*x + 0*x^2)*(0.1555338)+(0.1586638)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.2370150 + 0*x + 0*x^2)*0.1555338 + 0.1586638}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.2370150 + 0*x + 0*x^2)*0.1555338 + 0.1586638}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.2370150 + 0*x + 0*x^2)*0.1555338 + 0.1586638}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.2370150 + 0*x + 0*x^2)*0.1555338 + 0.1586638}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.2370150 + (0.0471861000+0.9040450)*x + (-0.0020322400-0.0741831)*x^2)*0.1555338 + 0.1586638}, size=3, xlim=c(0,8), colour='#EC1804')

# print(richnessPlot) #export at 1200x1000






###evenness
chainsEvennessIntercept <- chainsCommunity[,585:874]%>%
  gather(key=parameter, value=value, B.1.3.1:B.289.3.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsEvennessIntercept)[1] <- 'parameter1'

chainsEvennessSlope <- chainsCommunity[,1741:2029]%>%
  gather(key=parameter, value=value, B.1.3.2:B.289.3.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsEvennessSlope)[1] <- 'parameter2'

chainsEvennessQuad <- chainsCommunity[,2897:3185]%>%
  gather(key=parameter, value=value, B.1.3.3:B.289.3.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-2*sd, upper=quad+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
  select(parameter, quad)
names(chainsEvennessQuad)[1] <- 'parameter3'


#merge together with experiment list
evenness <- cbind(chainsEvennessIntercept, chainsEvennessSlope, chainsEvennessQuad)%>%
  select(-parameter2, -parameter3)%>%
  separate(parameter1, into=c('head', 'exp', 'type', 'variable'), sep='\\.', remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

evenness2 <- evenness[order(evenness$exp2),]
evenness3 <- cbind(expInfo2, evenness2)%>%
  select(-exp2, -exp)

evenness4 <- left_join(evenness3, expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  #get estimates at various time points
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.0739497+0.07095724,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.0739497+0.07095724,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.0739497+0.07095724)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(evenness4,'evenness_equations.csv', row.names=F)


#main figure
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(0,0.6))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.698156*x + -0.06299575*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.524705 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4598315 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5996665 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5993795 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.46968 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4559865 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.468053 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4893705 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6147075 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.506871 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5303835 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6354155 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.459786*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.796375*x + -0.1457845*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5191725 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.6177015*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.6417765 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.573123 + 0.425522*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.581519 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.583074 + 0.569577*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 2.339715*x + -0.2004035*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.577156 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.579755 + 0.9482115*x + -0.0625094*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(-0.593993 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.5607065 + 0.679758*x + -0.04949015*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-0.632632 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.530004 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3761605*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.423552 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.455183 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4177945 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4414265 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.48108 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5107825 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.49134 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.508146 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4616925*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.074075 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.9963765 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.9827085 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.794683*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.987199*x + -0.0684385*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 1.256075*x + -0.0755236*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4720895*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.365811*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.530052*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5263965 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5168835 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.537706 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.546505 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2676605*x + -0.0190846*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.520829 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.494171*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0739497)+(0.07095724)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.2847305 + 0*x + 0*x^2)*0.0739497 + 0.07095724}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.2847305 + 0*x + 0*x^2)*0.0739497 + 0.07095724}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.2847305 + 0*x + 0*x^2)*0.0739497 + 0.07095724}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.2847305 + 0*x + 0*x^2)*0.0739497 + 0.07095724}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.2847305 + (0.0751436000+0.8376550)*x + (-0.0008198390-0.0622513)*x^2)*0.0739497 + 0.07095724}, size=3, xlim=c(0,8), colour='#EC1804')

# print(evennessPlot) #export at 1200x1000


#print all plots together
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000




# ###density plots of all raw data
# meanDensity <- ggplot(data=rawData, aes(x=mean_change)) +
#   geom_density() +
#   xlab('Mean Change') +
#   ylab('Density') +
#   xlim(0,1) +
#   ylim(0,6.3)
# dispersionDensity <- ggplot(data=rawData, aes(x=dispersion_change)) +
#   geom_density() +
#   xlab('Dispersion Change') +
#   ylab('') +
#   geom_vline(xintercept=0, lty=2) +
#   xlim(-1,1.25) +
#   ylim(0,6.3)
# richnessDensity <- ggplot(data=rawData, aes(x=S_PC)) +
#   geom_density() +
#   xlab('Proportion Richness Change') +
#   ylab('Density') +
#   geom_vline(xintercept=0, lty=2) +
#   xlim(-1,1.25) +
#   ylim(0,6.3)
# evennessDensity <- ggplot(data=rawData, aes(x=SimpEven_change)) +
#   geom_density() +
#   xlab('Evenness Change') +
#   ylab('') +
#   geom_vline(xintercept=0, lty=2) +
#   xlim(-1,1.25) +
#   ylim(0,6.3)
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanDensity, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(dispersionDensity, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(richnessDensity, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(evennessDensity, vp=viewport(layout.pos.row=2, layout.pos.col=2))
# #export at 1200 x 1000
# 
# 
# ###density plots of final year estimates
# meanFinalDensity <- ggplot(data=mean4, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Mean Change') +
#   ylab('Density') +
#   xlim(0,1) +
#   ylim(0,10)
# dispersionFinalDensity <- ggplot(data=dispersion4, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Dispersion Change') +
#   ylab('') +
#   geom_vline(xintercept=0, lty=2) +
#   xlim(-1,1.25) +
#   ylim(0,10)
# richnessFinalDensity <- ggplot(data=richness4, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Proportion Richness Change') +
#   ylab('Density') +
#   geom_vline(xintercept=0, lty=2) +
#   xlim(-1,1.25) +
#   ylim(0,10)
# evennessFinalDensity <- ggplot(data=evenness4, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Evenness Change') +
#   ylab('') +
#   geom_vline(xintercept=0, lty=2) +
#   xlim(-1,1.25) +
#   ylim(0,10)
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanFinalDensity, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(dispersionFinalDensity, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(richnessFinalDensity, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(evennessFinalDensity, vp=viewport(layout.pos.row=2, layout.pos.col=2))
# #export at 1200 x 1000







###by resource mani
#still need to calculate the proportion of chains where x resource response was greater than y resource response

#mean change
meanResourceDrought <- mean4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
meanResource <- mean4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(meanResourceDrought)

#dispersion change
dispersionResourceDrought <- dispersion4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
dispersionResource <- dispersion4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(dispersionResourceDrought)

#richness change
richnessResourceDrought <- richness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
richnessResource <- richness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(richnessResourceDrought)

#evenness change
evennessResourceDrought <- evenness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
evennessResource <- evenness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr9, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(evennessResourceDrought)


#by resource at year 9
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=meanResource, variable='yr9', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(0, 0.41)) +
  xlab('')+
  annotate('text', x=0.5, y=0.40, label='(a)', size=10, hjust='left')
dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=dispersionResource, variable='yr9', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.08, 0.08, 0.04), name='Change in Dispersion') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  # coord_cartesian(ylim=c(-0.08, 0.08)) +
  xlab('')+
  annotate('text', x=0.5, y=0.038, label='(b)', size=10, hjust='left')
richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=richnessResource, variable='yr9', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  # coord_cartesian(ylim=c(-0.11, 0.1)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.1, label='(c)', size=10, hjust='left')
evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=evennessResource, variable='yr9', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.05, 0.05, 0.02), name='Change in Evenness') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  # coord_cartesian(ylim=c(0, 0.045)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.05, label='(d)', size=10, hjust='left')

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
  ylim(-0.5, 1.15) +
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
  ylim(-0.2, 0.7) +
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
  ylim(-0.055, 0.02) +
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
  ylim(-0.4, 1.4) +
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
  ylim(-0.4, 0.5) +
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
  ylim(-0.05, 0.05) +
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
  ylim(-0.4, 0.8) +
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
  ylim(-0.2, 1.2) +
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
  ylim(-0.12, 0.04) +
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
  ylim(-0.5, 0.5) +
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
  ylim(-0.3, 1.2) +
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
  ylim(-0.1, 0.04) +
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
  scale_y_continuous(limits=c(-0.6, 0.2), breaks=seq(-1, 0.2, 0.2)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Dispersion\nChange') +
  annotate('text', x=3.45, y=-0.6, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='richness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.4, 0.2), breaks=seq(-0.4, 0.4, 0.2)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Proportion\nRichness Change') +
  annotate('text', x=3.45, y=-0.4, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.5, 0.25), breaks=seq(-0.4, 0.3, 0.2)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.5, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
#export at 2400x500

# 
# 
# 
# 
# # ###look for patterns of spp appearance/disappearance -- no clear patterns, probably because just the few CDR examples that are long term enough to see the pattern
# relAbund <- read.csv('SpeciesRelativeAbundance_April2016.csv')%>%
#   select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
#   mutate(exp_trt=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
#   #get rid of duplicate species within a plot and year in the dataset; once we contact the dataowners, this step will no longer be needed
#   group_by(exp_trt, site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species)%>%
#   summarise(relcov=mean(relcov))%>%
#   filter(exp_trt!='NIN::herbdiv::0::5F' & site_code!='GVN')
# 
# expinfo<-read.csv('ExperimentInformation_Mar2016.csv')%>%
#   mutate(exp_trt=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
#   select(exp_trt, plot_mani, calendar_year)
# 
# relAbundYear<-merge(relAbund, expinfo, by=c("exp_trt","calendar_year"), all=F)
# 
# #make a new dataframe with just the label
# exp_trt=relAbundYear%>%
#   select(exp_trt)%>%
#   unique()
# 
# #make a new dataframe to collect the turnover metrics
# turnoverAll=data.frame(row.names=1)
# 
# for(i in 1:length(relAbundYear$exp_trt)) {
# 
#   #creates a dataset for each unique year, trt, exp combo
#   subset=relAbundYear[relAbundYear$exp_trt==as.character(exp_trt$exp_trt[i]),]%>%
#     select(exp_trt, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
#     #get just first and last year of study
#     filter(calendar_year==min(calendar_year)|calendar_year==max(calendar_year))
# 
#   #need this to keep track of plot mani
#   labels=subset%>%
#     select(exp_trt, plot_mani, calendar_year)%>%
#     unique()
# 
#   #calculate disappearance
#   disappearance=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', replicate.var=NA, metric='disappearance')%>%
#     group_by(calendar_year)%>%
#     summarise(disappearance=mean(disappearance))
# 
#   #calculate appearance
#   appearance=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', replicate.var=NA, metric='appearance')%>%
#     group_by(calendar_year)%>%
#     summarise(appearance=mean(appearance))
#   
#   #calculate turnover
#   total=turnover(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', replicate.var=NA, metric='total')%>%
#     group_by(calendar_year)%>%
#     summarise(turnover=mean(total))
# 
#   #merging back with labels to get back plot_mani
#   turnover=labels%>%
#     left_join(disappearance, by='calendar_year')%>%
#     left_join(appearance, by='calendar_year')%>%
#     left_join(total, by='calendar_year')%>%
#     filter(calendar_year==max(calendar_year))%>%
#     select(exp_trt, plot_mani, appearance, disappearance, turnover)
# 
#   #pasting variables into the dataframe made for this analysis
#   turnoverAll=rbind(turnover, turnoverAll)
# }
# 
# turnoverCtl <- turnoverAll%>%
#   filter(plot_mani==0)%>%
#   separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)%>%
#   select(site_code, project_name, community_type, appearance, disappearance, turnover)
# names(turnoverCtl)[names(turnoverCtl)=='appearance'] <- 'appearance_ctl'
# names(turnoverCtl)[names(turnoverCtl)=='disappearance'] <- 'disappearance_ctl'
# names(turnoverCtl)[names(turnoverCtl)=='turnover'] <- 'turnover_ctl'
# 
# turnoverDiff <- turnoverAll%>%
#   mutate(trt=ifelse(plot_mani==0, 'ctl', 'trt'))%>%
#   separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)%>%
#   filter(trt!='ctl')%>%
#   left_join(turnoverCtl, by=c('site_code', 'project_name', 'community_type'))%>%
#   mutate(appearance_diff=appearance-appearance_ctl, disappearance_diff=disappearance-disappearance_ctl, turnover_diff=turnover-turnover_ctl)
# 
# plot(turnoverDiff$plot_mani, turnoverDiff$appearance_diff)
# plot(turnoverDiff$plot_mani, turnoverDiff$disappearance_diff)
# plot(turnoverDiff$plot_mani, turnoverDiff$turnover_diff)
# plot(turnoverDiff$plot_mani, turnoverDiff$appearance)
# plot(turnoverDiff$plot_mani, turnoverDiff$disappearance)
# plot(turnoverDiff$plot_mani, turnoverDiff$turnover)
# 
# summary(glm(turnover~as.factor(plot_mani), data=turnoverDiff))
# lsmeans(glm(turnover~as.factor(plot_mani), data=turnoverDiff), 'plot_mani')
# 
# ggplot(data=turnoverDiff, aes(x=as.factor(plot_mani), y=turnover)) +
#   geom_boxplot() +
#   xlab('Number of Factors Manipulated') +
#   ylab('Species Turnover')
# 
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
# 
# 
# 
# 
# # ###look at spp comp of five factor manipulations to find patterns of immigration or loss of dominant spp
# # relAbundFive <- relAbundYear
# # 
# # #make a new dataframe with just the label
# # expTrtYear=relAbundFive%>%
# #   select(exp_trt)%>%
# #   unique()
# # 
# # #make a new dataframe to collect the turnover metrics
# # relAbundFiveYear=data.frame(row.names=1)
# # 
# # for(i in 1:length(expTrtYear$exp_trt)) {
# # 
# #   #creates a dataset for each unique year, trt, exp combo
# #   subset=relAbundFive[relAbundFive$exp_trt==as.character(expTrtYear$exp_trt[i]),]%>%
# #     select(exp_trt, calendar_year, treatment, plot_mani, genus_species, relcov, plot_id)%>%
# #     group_by(exp_trt, calendar_year, treatment, plot_mani, genus_species)%>%
# #     summarise(relcov=mean(relcov))%>%
# #     ungroup()%>%
# #     #get just first and last year of study
# #     filter(calendar_year==min(calendar_year)|calendar_year==max(calendar_year))%>%
# #     mutate(time=ifelse(calendar_year==min(calendar_year), 'first', 'last'))%>%
# #     select(-calendar_year, -treatment)%>%
# #     spread(key=time, value=relcov, fill=0)
# # 
# #   #pasting variables into the dataframe made for this analysis
# #   relAbundFiveYear=rbind(subset, relAbundFiveYear)
# # }
# # 
# # relAbundFiveYear <- relAbundFiveYear%>%
# #   separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)
# # 
# # relAbundFiveYearRich <- richness4%>%
# #   left_join(relAbundFiveYear, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
# #   select(site_code, project_name, community_type, treatment, plot_mani, genus_species, first, last, experiment_length, intercept, slope, quad, nutrients, water, carbon, precip, alt_length, yr9)%>%
# #   filter(plot_mani>4)
# # 
# # ggplot(data=relAbundFiveYearRich, aes(x=first, y=last, colour=quad)) +
# #   geom_point() +
# #   scale_colour_gradientn(colours=rainbow(4))
# 
# 
# 
# 
# 
# #look at number replicates for dispersion results -- doesn't make a difference
# reps <- relAbund%>%
#   group_by(site_code, project_name, community_type, treatment, calendar_year, plot_id)%>%
#   summarise(mean=mean(relcov))%>%
#   ungroup()%>%
#   group_by(site_code, project_name, community_type, treatment, calendar_year)%>%
#   summarise(rep_num=n())%>%
#   ungroup()%>%
#   group_by(site_code, project_name, community_type, treatment)%>%
#   summarise(rep_num=mean(rep_num))
# 
# dispersionReps <- dispersion4%>%
#   left_join(reps, by=c('site_code', 'project_name', 'community_type', 'treatment'), all=F)%>%
#   #filter out lovegrass mistake trt code until fixed in main dataset
#   filter(rep_num>1)
# 
# # ggplot(data=dispersionReps, aes(x=rep_num, y=intercept, colour=plot_mani)) +
# #   geom_point()
# # ggplot(data=dispersionReps, aes(x=rep_num, y=slope, colour=plot_mani)) +
# #   geom_point()
# # ggplot(data=dispersionReps, aes(x=rep_num, y=quad, colour=plot_mani)) +
# #   geom_point()
# # ggplot(data=dispersionReps, aes(x=rep_num, y=yr10, colour=plot_mani)) +
# #   geom_point()
# 
# ggplot(data=dispersionReps, aes(x=rep_num, y=yr9)) +
#   geom_point() +
#   xlab('Number of Relicates') +
#   ylab('Dispersion Change') +
#   scale_x_continuous(breaks=seq(0,50,10)) +
#   coord_cartesian(xlim=c(0,45))
# #export at 900x900
#   
# 
# 
# 
# 
# ###look at five factor manipulations for mean change
# # #just for the four experiments with five factors, compare to their four factor treatments
# # meanFive <- mean4%>%
# #   filter(treatment=='1_y_n'|treatment=='8_y_n'|treatment=='1_f_u_n'|treatment=='8_f_u_n'|treatment=='2F'|treatment=='3F'|treatment=='4F'|treatment=='ghn'|treatment=='gsn'|treatment=='ncn'|treatment=='nhn'|treatment=='nsn')
# # 
# # cdr1APlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='A'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
# #   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
# #                      labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(a) CDR e001 A', size=10, hjust='left') +
# #   theme(legend.position='none')
# # cdr1BPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='B'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
# #                    labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(b) CDR e001 A', size=10, hjust='left') +
# #   theme(legend.position='none')
# # cdr1CPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='C'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
# #                    labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(c) CDR e001 A', size=10, hjust='left') +
# #   theme(legend.position='none')
# # cdr1DPlot <- ggplot(data=subset(meanFive, project_name=='e001'&community_type=='A'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('1_y_n', '8_y_n'),
# #                    labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(d) CDR e001 D', size=10, hjust='left') +
# #   theme(legend.position='none')
# # ninPlot <- ggplot(data=subset(meanFive, site_code=='NIN'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('2F', '3F', '4F'),
# #                    labels=c('4a', '4b', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(e) NIN herbdiv', size=10, hjust='left') +
# #   theme(legend.position='none')
# # cdr2APlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='A'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
# #   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
# #                    labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(f) CDR e002 A', size=10, hjust='left') +
# #   theme(legend.position='none')
# # cdr2BPlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='B'), aes(x=treatment, y=yr10, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
# #                    labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(g) CDR e002 B', size=10, hjust='left') +
# #   theme(legend.position='none')
# # cdr2CPlot <- ggplot(data=subset(meanFive, project_name=='e002'&community_type=='C'), aes(x=treatment, y=final_year_estimate, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('1_f_u_n', '8_f_u_n'),
# #                    labels=c('4', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(h) CDR e002 C', size=10, hjust='left') +
# #   theme(legend.position='none')
# # traPlot <- ggplot(data=subset(meanFive, site_code=='TRA'), aes(x=treatment, y=final_year_estimate, fill=treatment)) +
# #   geom_bar(stat="identity", colour='black') +
# #   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
# #   scale_x_discrete(limits=c('ghn', 'gsn', 'ncn', 'nhn', 'nsn'),
# #                    labels=c('4', '4', '4', '5', '5')) +
# #   coord_cartesian(ylim=c(0,1)) +
# #   scale_fill_manual(values=c('white', 'white', 'white', 'black', 'black')) +
# #   xlab('') +
# #   annotate('text', x=0.5, y=1, label='(i) TRA lovegrass', size=10, hjust='left') +
# #   theme(legend.position='none')
# # 
# # pushViewport(viewport(layout=grid.layout(2,5)))
# # print(cdr1APlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# # print(cdr1BPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# # print(cdr1CPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
# # print(cdr1DPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
# # print(ninPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 5))
# # print(cdr2APlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# # print(cdr2BPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# # print(cdr2CPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
# # print(traPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
# # #export at 2400x1200
# 
# 
# #compare any four factor without N to five factor with N
# expRawMean <- expRaw%>%
#   group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
#   summarise(n=mean(n), herb_removal=mean(herb_removal), plant_mani=mean(plant_mani))
# 
# meanCompare <- mean4%>%
#   left_join(expRawMean, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
#   mutate(n_mani=ifelse(n>0, 1, 0))
# 
# 
# #plot without N at four factors, with N at five factors
# compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
#   geom_bar(stat="identity", colour='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
#   scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'grey')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
#   theme(legend.position='none')
# compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&herb_removal>0), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
#   geom_bar(stat="identity", colour='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'grey')) +
#   xlab('Number of Factors Manipulated') +
#   annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
#   theme(legend.position='none')
# comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plant_mani>0), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
#   geom_bar(stat="identity", colour='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
#   scale_x_discrete(labels=c('4 factor\n+manip.', '5 factor\n+manip.')) +
#   coord_cartesian(ylim=c(0,1)) +
#   scale_fill_manual(values=c('white', 'grey')) +
#   xlab('') +
#   annotate('text', x=0.5, y=1, label='(c) Plant Manipulation Comparison', size=10, hjust='left') +
#   theme(legend.position='none')
# 
# pushViewport(viewport(layout=grid.layout(1,3)))
# print(compareNPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(comparePlantPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
# print(compareHerbPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# #export at 2400x1200
# 
# 
# 
# 
