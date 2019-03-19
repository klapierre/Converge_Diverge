library(ggplot2)
library(grid)
library(mgcv)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

setwd('C:\\Users\\Kim\\Desktop\\bayesian output\\anpp site subset')

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
setwd('C:\\Users\\Kim\\Desktop\\bayesian output')

#experiment information
expRaw <- read.csv('ExperimentInformation_Mar2016.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), precip=mean(precip))

rawData <- read.csv('ForBayesianAnalysis_CommWithANPP_July2016.csv')

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
chains1 <- read.csv('diversity_MVANPP_rdisp_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('diversity_MVANPP_rdisp_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('diversity_MVANPP_rdisp_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('diversity_MVANPP_rdisp_3.csv', comment.char='#')
chains4 <- chains4[-1:-5000,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)

#get values for overall (mean) lines across levels of plot mani
#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4)%>%
  gather(key=parameter, value=value, U_int.2.1:mu_quad.4)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, keep=ifelse(diff==-2, 0, median ))

###mean change
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsMeanIntercept <- chainsCommunity[,7:183]%>%
  gather(key=parameter, value=value, B.1.1.1:B.177.1.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsMeanIntercept)[1] <- 'parameter1'

chainsMeanSlope <- chainsCommunity[,715:891]%>%
  gather(key=parameter, value=value, B.1.1.2:B.177.1.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsMeanSlope)[1] <- 'parameter2'

chainsMeanQuad <- chainsCommunity[,1423:1599]%>%
  gather(key=parameter, value=value, B.1.1.3:B.177.1.3)%>%
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
  #get estimates at various time points
  mutate(yr10=(intercept + 10*slope + (10^2)*quad)*0.191162+0.3322108,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.191162+0.3322108,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.191162+0.3322108)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(mean4,'mean_equations.csv', row.names=F)



#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  ylim(-1,2) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.24691*x + -0.012965*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.509469 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0.5978675 + 0.196978*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.8872985 + 0.227538*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.002485 + 0.277576*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.8176425 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7232765 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8998685 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4157465 + 0.163817*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.5766025 + 0.16471*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.4735015 + 0.162621*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.3830325 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2346875*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.744456 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6787955 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.944192 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6494645 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.55009 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1755535*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.170387*x + -0.01999125*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5778505 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4958605 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6033275 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.523705 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.421004 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.426122 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.192615*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.515272 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.487531 + 0.1951095*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6640915 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.187766*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.258705*x + -0.02005925*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.709477 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.764141 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7448235 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(-1.27029 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0.451995 + 0.285277*x + -0.00934287*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.899032 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.3935695 + 0.45654*x + -0.01381745*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5473395 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0.5545965 + 0.178807*x + -0.00622573*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + -0.00636409*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(1.09447 + 0.262483*x + -0.01405435*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7106665 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2571585*x + -0.01253035*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.4067645 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0.8039855 + 0.4811675*x + -0.0225637*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.576431 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4100905*x + -0.015939*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.739169 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6136435 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7928375 + 0.21809*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9995195 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.432499 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.455843 + 0.2826325*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8152595 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.563572 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(1.620955 + 0.1303165*x + -0.004276105*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5174705 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.611502 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.786625 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.1578145*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.1526275*x + -0.00588756*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.1971515*x + -0.00917203*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.129583*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.629245 + 0.11365*x + -0.00722325*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6067495 + 0.1522635*x + -0.00972466*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.117668*x + -0.00634147*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(-0.776657 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.79809 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.693897 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.825771 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.999503 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0518 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.008835 + 0.126748*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8861555 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.746136 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5777045 + 0.180061*x + -0.00899127*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.16765*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.755297 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.952752 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.873431 + 0.1254205*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.984447 + 0.166349*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.502106 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8017685 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.762147 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7426495 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2293495*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.4550925 + 0.2242875*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.919753 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6295045 + 0.2230595*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.657429 + 0.185564*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.611999 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.546951 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.190099*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.549558 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.504108 + 0.186006*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.506621 + 0.232495*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.719328 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7662025 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-1.080265 + 0.1535975*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.0925 + 0.1889185*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8141845 + 0.3674945*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0827 + 0.167018*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.32201 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.16079 + 0.2203485*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.176995 + 0.1469045*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.44606 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.743094 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.620448 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6331315 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7602405 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.745985 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.745556 + 0*x + -0.01287355*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.34278 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.230725 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9909945 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-1.558435 + 0.1789755*x + -0.009611065*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-0.441973 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9653275 + 0.1745495*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9552315 + 0.19189*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.897426 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9270975 + 0.1786065*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.719565 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8549055 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.2783335*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9779835 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.559003 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.191162)+(0.3322108)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.437075 + 0.103359*x + 0*x^2)*0.191162 + 0.3322108}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.437075 + 0.103359*x + 0*x^2)*0.191162 + 0.3322108}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.437075 + 0.103359*x + 0*x^2)*0.191162 + 0.3322108}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.437075 + 0.103359*x + 0*x^2)*0.191162 + 0.3322108}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.437075 + 0.7151795) + (0.103359 + 0.25281)*x + (-0.00582201 + -0.009176025)*x^2)*0.191162 + 0.3322108}, size=3, xlim=c(0,22), colour='#EC1804')

print(meanPlot) #export at 1200x1000