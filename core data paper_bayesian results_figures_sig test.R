library(ggplot2)
library(grid)
library(mgcv)
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

rawData <- read.csv('ForBayesianAnalysis_March2016b.csv')

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
chainsCommunity <- read.csv('fullChains.csv')

#get values for overall (mean) lines across levels of plot mani
#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4)%>%
  gather(key=parameter, value=value, U_int.2.1:mu_quad.4)%>%
  group_by(parameter)%>%
  summarise(median=median(value))

###mean change
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsMeanIntercept <- chainsCommunity[,7:296]%>%
  gather(key=parameter, value=value, B.1.1.1:B.290.1.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsMeanIntercept)[1] <- 'parameter1'

chainsMeanSlope <- chainsCommunity[,1167:1456]%>%
  gather(key=parameter, value=value, B.1.1.2:B.290.1.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsMeanSlope)[1] <- 'parameter2'

chainsMeanQuad <- chainsCommunity[,2327:2616]%>%
  gather(key=parameter, value=value, B.1.1.3:B.290.1.3)%>%
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
  mutate(yr10=(intercept + 10*slope + (10^2)*quad)*0.1701297+0.3140121,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.1701297+0.3140121,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.1701297+0.3140121)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,',
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
  stat_function(fun=function(x){(-0.4883495 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.631301 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5098295 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2023475*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2527795*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5530665 + 0.3011225*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.331224*x + -0.01298095*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.7045975 + 0.383736*x + -0.015702*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.27143 + 0.318652*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.12772 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.02848 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.00456 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.42911 + 0.2339005*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.577968 + 0.210579*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.4461505 + 0.206251*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.248163*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.582217 + 0.266819*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2203325*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8106405 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.683987 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.02411 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7965215 + 0.202106*x + -0.0113552*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-0.609518 + 0.260602*x + -0.0151973*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8540235 + 0.1386425*x + -0.00978064*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8374905 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-1.51362 + 0.184485*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.63217 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.59007 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.63067 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.50476 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.662305 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6354955 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5043855 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1745365*x + -0.01733725*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5628235 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4731355 + 0.1710675*x + -0.017348*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6129605 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4940195 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1837785*x + -0.01724645*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1867935*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6138175 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1762135*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2457725*x + -0.0170403*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.669725 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.671122 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.707299 + 0.1289235*x + -0.00888552*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(-1.30051 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0.589338 + 0.3214405*x + -0.01045525*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.8919455 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.445834 + 0.530661*x + -0.01611715*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.4873165 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0.6912785 + 0.2038135*x + -0.00703984*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.1335225*x + -0.008427285*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(1.248975 + 0.3036995*x + -0.01604555*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5986995 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2480395*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(1.026485 + 0.543711*x + -0.0258208*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5274785 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4682845*x + -0.0170641*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7563585 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.57131 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7582515 + 0.2581355*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.05297 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-1.348815 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.24015 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.19757*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4165555 + 0.176853*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4400815 + 0.3280515*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.844727 + 0.2267175*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.269437*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.222553*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7693165 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.544599 + 0.182794*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.63166 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.80487 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.85503 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.902496 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.193055 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.014385 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.103805 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.127005 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9711095 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.590733 + 0.189764*x + -0.01111425*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.400373 + 0.137916*x + -0.009604315*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4380145 + 0.173627*x + -0.0108367*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4158195 + 0.127494*x + -0.008657875*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.149924*x + -0.0100193*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.455317 + 0.2107825*x + -0.0133669*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5207225 + 0.193823*x + -0.01140495*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.198613*x + -0.0127062*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.521741 + 0.196093*x + -0.01131145*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.516351 + 0.160993*x + -0.01019485*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.416823 + 0.180059*x + -0.01064965*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.646893 + 0.1334485*x + -0.0083625*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5936595 + 0.111935*x + -0.0080649*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.629567 + 0.132426*x + -0.009647625*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.429283 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(1.956205 + 0.142431*x + -0.004659135*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8367625 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.692884 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9020235 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.521387 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.126315*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.194454*x + -0.006083915*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.20106*x + -0.00794801*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.245699*x + -0.01145405*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1260425*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.814752 + 0.13187*x + -0.008340265*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.773781 + 0.195075*x + -0.0122432*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.125592*x + -0.006563*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(-0.705988 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7657325 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6265235 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8095375 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9669255 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.07371 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.05278 + 0.1521105*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.863773 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.75 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5707205 + 0.2137465*x + -0.0106657*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.261015*x + -0.0129277*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.815515 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00993 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.897059 + 0.165041*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.002265 + 0.193276*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.639673 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4271265 + 0.1881695*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.184442*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1861185*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.270301*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3308325*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.313874*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.215788*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.226289*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4703485 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8601045 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.799565 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7766935 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.8730395 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.5780865 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.503424 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.415608 + 0.2575895*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2517675*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5701995 + 0.231955*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7925675 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.807409 + 0.1946905*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.830252 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.647873 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.947481 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6842465 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.226134*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.939811 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6155325 + 0.2556765*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.630523 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.617254 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5125785 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2205565*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.448908 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.69893 + 0.262698*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.459635 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.619206 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.727828 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3191835*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.2791165*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2810145*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3150845*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(-0.838739 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.15033 + 0.1968865*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7507745 + 0.386003*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.11676 + 0.178316*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.298205 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.12568 + 0.1978265*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.11508 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9910815 + 0.1844125*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-1.001215 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9765795 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.013995 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.912906 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8685015 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4604675 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.456075 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.67413 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.526694 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.541302 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.747209 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7339715 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4589845 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.579963 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.533806 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.684806 + 0*x + -0.0146405*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6100035 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.390845 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.329735 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.91243 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-1.649315 + 0.2199335*x + -0.01214475*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-0.53718 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.076905 + 0.234262*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.083185 + 0.246818*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.0608 + 0.1684815*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.062665 + 0.2154275*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6313155 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.847431 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.753947 + 0.220392*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8711005 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.549339 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2732865*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.921255 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5786435 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4673615 + 0.182799*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.740998 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.511347 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.478088 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.199298*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8887615 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.622662 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.424854 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0.626957 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.5441655 + 0.18185) + (0.132175)*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.5441655 + 0.731963) + (0.132175 + 0.212783)*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#EC1804')

# print(meanPlot) #export at 1200x1000




###dispersion
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsDispersionIntercept <- chainsCommunity[,297:586]%>%
  gather(key=parameter, value=value, B.1.2.1:B.290.2.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsDispersionIntercept)[1] <- 'parameter1'

chainsDispersionSlope <- chainsCommunity[,1457:1746]%>%
  gather(key=parameter, value=value, B.1.2.2:B.290.2.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsDispersionSlope)[1] <- 'parameter2'

chainsDispersionQuad <- chainsCommunity[,2617:2906]%>%
  gather(key=parameter, value=value, B.1.2.3:B.290.2.3)%>%
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
  #get estimates at various time points
  mutate(yr10=(intercept + 10*slope + (10^2)*quad)*0.09064568-0.00235573,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.09064568-0.00235573,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.09064568-0.00235573)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, 'grey', ifelse(plot_mani==2, 'grey', ifelse(plot_mani==3, 'grey', ifelse(plot_mani==4, 'grey', 'grey')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(dispersion4,'dispersion_equations.csv', row.names=F)



#main figure
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.35,0.35))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.35, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.733904 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.01589495*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.01763095*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02013605*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.31945*x + 0.0483592*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.046109*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.6929055 + -0.2895485*x + 0.04496775*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0444857*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.04136 + 0*x + 0.04483635*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3781695*x + 0.0518324*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.337569*x + 0.0494296*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3268365*x + 0.0490314*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.713692 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7435525 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7272455 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.8585995 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(-2.2985 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(1.09329 + -0.1727785*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.006400485*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0.759779 + 0.197894*x + -0.01064735*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,20), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,20), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.51141 + -0.2921915*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.390211*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.491978*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.220695 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.8419925 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-1.98304 + 0.502853*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.93097 + 0.4285995*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.50593 + 0.4089585*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.123815 + 0.4019995*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.44016 + 0.508324*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.9833 + 0.4921935*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.814166 + 0.4588745*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.1773 + 0.443874*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3100105*x + 0.01105735*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.01207 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.18472*x + 0.008728115*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.245089*x + 0.0117654*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.302262*x + 0.0143793*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1971325*x + 0.00869624*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.197901*x + 0.00933956*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.793563 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1935165*x + 0.008523125*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.21274 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1622775*x + 0.01091985*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.13541 + 0.3290465*x + -0.029531*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.19515 + 0*x + -0.02142695*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.02260125*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0184907*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.897932 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7039625 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.257568*x + 0.03824745*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.849515 + -0.3758145*x + 0.03894875*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3386565*x + 0.0404048*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0576143*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0600158*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.276059*x + 0.0639807*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4161715*x + 0.0688968*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0353412*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2832695*x + 0.0383121*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0352133*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0370573*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2234215*x + 0.03708675*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03639645*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03340775*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03352205*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03779605*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03298385*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.000405 + 0*x + 0.03249165*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.893366 + 0*x + 0.0355578*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9058435 + 0*x + 0.03444065*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.958849 + 0*x + 0.03662495*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.307535 + 0*x + 0.0321784*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.145325 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.67039 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.64099 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7713325 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.839616 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.14491 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6590865 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6522305 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3031565*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9592215 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.733284 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.03674 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7710045 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.712423 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.12243 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.748319 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.36553 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8500545 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.484455 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.422335 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6893565 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  
#estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0 - 0}, size=3, xlim=c(0,23), colour='black')

# print(dispersionPlot) #export at 1200x1000







###richness
chainsRichnessIntercept <- chainsCommunity[,877:1166]%>%
  gather(key=parameter, value=value, B.1.4.1:B.290.4.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsRichnessIntercept)[1] <- 'parameter1'

chainsRichnessSlope <- chainsCommunity[,2037:2326]%>%
  gather(key=parameter, value=value, B.1.4.2:B.290.4.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsRichnessSlope)[1] <- 'parameter2'

chainsRichnessQuad <- chainsCommunity[,3197:3486]%>%
  gather(key=parameter, value=value, B.1.4.3:B.290.4.3)%>%
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
  #get estimates at various time points
  mutate(yr10=(intercept + 10*slope + (10^2)*quad)*0.2287037-0.07758351,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.2287037-0.07758351,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.2287037-0.07758351)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(richness4,'richness_equations.csv', row.names=F)



#main figure
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-1,0.65))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change') +
  annotate('text', x=0, y=0.65, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6738265 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0.6948925 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2268595*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.612081 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.297744*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3895715*x + 0.0168034*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2827385*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5514035 + -0.3887175*x + 0.0484527*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0347035*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.03711085*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.03368555*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.336831*x + 0.04628665*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.319147*x + 0.0439466*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.3267065*x + 0.04505255*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6272095 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5268475 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5362985 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.239153*x + 0.014227*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5491785 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.9366445 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.5430035 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.6874735 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.763295 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8738605 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.5635785 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.295472*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.612026 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.238211*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.22296*x + 0.0110055*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.096995 + -0.36493*x + 0.01945115*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.879221 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.3347215*x + 0.0105112*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.4786015 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.8939105 + -0.3195775*x + 0.0142767*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.912655 + -0.2791495*x + 0.0117167*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.21914*x + 0.0152463*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.11077 + -0.5381785*x + 0.0308385*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.576031 + -0.531151*x + 0.0279142*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.996043 + -0.424857*x + 0.02612505*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0.7338375 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6194855 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2148725*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.19861 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.72591 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.393202*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.595457 + -0.3009115*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3033635*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.61224 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4030275*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0.876623 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.893063 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.523034 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7323775 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0109624*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.645894 + 0*x + -0.01144835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.502073 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0100508*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.487423 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4596655 + 0*x + -0.01185675*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.01086365*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.583547 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7541255 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.6630925 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6900255 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5400705 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1726115*x + 0.007361315*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.264*x + 0.02425245*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2706855*x + 0.01345105*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.195684*x + 0.01244695*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.168131*x + 0.01204785*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.390939*x + 0.0255996*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.268837*x + 0.0175611*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1754695*x + 0.01133415*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8523075 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7341665 + 0*x + -0.01509915*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.753422 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.256078*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.528627 + 0*x + -0.01789425*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.169129*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2675145*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.969191 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.03897 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8844005 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.04319 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7314515 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.9770775 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.6419825 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.032273*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0332422*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.03567035*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6882725 + -0.284502*x + 0.0367779*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.2655285*x + 0.0350543*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.03119915*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.02606895*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2727295*x + 0.04042885*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0.584382 + 0*x + 0.03140805*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.0313835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1929605*x + 0.0314701*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.02833215*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0317326*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.03795895*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.281549*x + 0.0391019*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0350707*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.270964*x + 0.04065335*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.03586595*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0.02956835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.03202525*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.0295487*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.056555 + 0*x + 0.0246246*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.624395 + 0.4319585*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.861475 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.616545 + 0.453692*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.3117205*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.351106*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4424085*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.699695 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2792245*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.04209 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.6727815 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.212907*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.2118165*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(1.023795 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5338895 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.822188 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.61744 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0.725389 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.7658975 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.5951095 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6706775 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.74666 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.561176 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2500135*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.243432*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.07507 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1648525*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0.850658 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.558802 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0.5918315 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6819755 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.598073 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.784434 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.759163 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8480005 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.688117 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.637682 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.846348 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.6611805 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7063845 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(1.20301 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.820268 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.522843 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5255085 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.976253 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(0.319741 + -0*x + 0*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(0.319741 + -0*x + 0*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(0.319741 + -0*x + 0*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(0.319741 + -0*x + 0*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((0.319741-0.694714) + (-0.0520446-0.3041105)*x + 0*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,22), colour='#EC1804')

# print(richnessPlot) #export at 1200x1000






###evenness
chainsEvennessIntercept <- chainsCommunity[,587:876]%>%
  gather(key=parameter, value=value, B.1.3.1:B.290.3.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsEvennessIntercept)[1] <- 'parameter1'

chainsEvennessSlope <- chainsCommunity[,1747:2036]%>%
  gather(key=parameter, value=value, B.1.3.2:B.290.3.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsEvennessSlope)[1] <- 'parameter2'

chainsEvennessQuad <- chainsCommunity[,2907:3196]%>%
  gather(key=parameter, value=value, B.1.3.3:B.290.3.3)%>%
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
  #get estimates at various time points
  mutate(yr10=(intercept + 10*slope + (10^2)*quad)*0.1034254+0.019179,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.1034254+0.019179,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.1034254+0.019179)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(evenness4,'evenness_equations.csv', row.names=F)


#main figure
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.35,0.55))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.55, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.2869265*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.648768 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6159045 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6039485 + -0.3012245*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2106535*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8076745 + -0.3056815*x + 0.01889215*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.250867*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.562218 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.659341 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.891012 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.667963 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8423015 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.020562*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0232075*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2131235*x + -0.0164914*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2085045*x + -0.01640775*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.486815*x + -0.0296503*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.342991*x + -0.01215835*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3129785*x + -0.0142065*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(1.40602 + 0.9065815*x + -0.06015*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.630119*x + -0.0334193*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4648995*x + -0.02317545*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.23086*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5317225 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.537868 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7911565 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.772084 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7362095 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.117303*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6513185 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6305825 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.01092495*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1356835*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1760245*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6155435 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.179416*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.359805 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.33328 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.30181 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.295233*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.523093*x + -0.0257636*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.574114*x + -0.0297145*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.862082*x + -0.04261665*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.3199705*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3661685*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3524815*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.334261*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.364569*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.412672*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4349475*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2257315*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.2377275*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2512495*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.203584*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4702325*x + -0.02763595*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.266361*x + 0.02586255*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3298555*x + 0.0300639*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.208304*x + 0.0210999*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3190975*x + 0.03012865*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6836215 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6673645 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.60223 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.538797 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.217861*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0.5550905 + 0.21163*x + -0.01412815*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2041825*x + -0.0215562*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.208523*x + -0.0203584*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2114845*x + -0.02102815*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.509345 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7058875 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.227075*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3098325*x + 0.02541585*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2699515*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2845285*x + 0.0234434*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5308165 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6127565 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5366045 + 0.5073425*x + -0.03066115*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2481755*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(-0.594625 + 0.471253*x + -0.0282798*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.02539825*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0240929*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.02535545*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.02475015*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.02915045*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0242973*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.167144 + 0*x + -0*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.167144 + 0*x + -0*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.167144 + 0*x + -0*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.167144 + 0*x + -0*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.167144) + (0.01842445+0.445834)*x + (-0.00178622-0.0236129)*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,22), colour='#EC1804')

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
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
meanResource <- mean4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(meanResourceDrought)

#dispersion change
dispersionResourceDrought <- dispersion4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
dispersionResource <- dispersion4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(dispersionResourceDrought)

#richness change
richnessResourceDrought <- richness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
richnessResource <- richness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(richnessResourceDrought)

#evenness change
evennessResourceDrought <- evenness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
evennessResource <- evenness4%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(evennessResourceDrought)


#by resource at final year
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=meanResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(0, 0.41)) +
  xlab('')+
  annotate('text', x=0.5, y=0.40, label='(a)', size=10, hjust='left')
dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=dispersionResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.04, 0.04)) +
  xlab('')+
  annotate('text', x=0.5, y=0.038, label='(b)', size=10, hjust='left')
richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=richnessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.11, 0.05)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.045, label='(c)', size=10, hjust='left')
evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=evennessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(0, 0.045)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.043, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600

# #10 year estimate - 10 years is the median experiment length
# meanResourcePlot10 <- ggplot(data=barGraphStats(data=meanResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(0, 0.5)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.5, label='(a)', size=10, hjust='left')
# dispersionResourcePlot10 <- ggplot(data=barGraphStats(data=dispersionResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.06, 0.10, 0.04), name='Dispersion Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.06, 0.10)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.10, label='(b)', size=10, hjust='left')
# richnessResourcePlot10 <- ggplot(data=barGraphStats(data=richnessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.5, 0.1, 0.1), name='Proportion Richness Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.2, 0.1)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.1, label='(c)', size=10, hjust='left')
# evennessResourcePlot10 <- ggplot(data=barGraphStats(data=evennessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.06, 0.06, 0.02), name='Evenness Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.06, 0.06)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.06, label='(d)', size=10, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlot10, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlot10, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlot10, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlot10, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600




      
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
  ylim(-0.5, 0.6) +
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
  ylim(-0.035, 0.05) +
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
  ylim(-0.5, 0.6) +
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
  ylim(-0.035, 0.05) +
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
  ylim(-0.5, 0.6) +
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
  ylim(-0.035, 0.05) +
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
  ylim(-0.5, 0.6) +
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
  ylim(-0.035, 0.05) +
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
  scale_y_continuous(limits=c(-0.8, 0.2), breaks=seq(-0.5, 0.5, 0.5)) +
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
  scale_y_continuous(limits=c(-0.2, 0.15), breaks=seq(-0.2, 0.2, 0.2)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Dispersion\nChange') +
  annotate('text', x=3.45, y=-0.2, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='richness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.12, 0.5), breaks=seq(-0.4, 0.4, 0.4)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Proportion\nRichness Change') +
  annotate('text', x=3.45, y=-0.12, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.31, 0.1), breaks=seq(-0.3, 0.3, 0.3)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.31, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
#export at 2400x500





# ###look for patterns of spp appearance/disappearance -- no clear patterns, probably because just the few CDR examples that are long term enough to see the pattern
relAbund <- read.csv('SpeciesRelativeAbundance_April2016.csv')%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_trt=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
  #get rid of duplicate species within a plot and year in the dataset; once we contact the dataowners, this step will no longer be needed
  group_by(exp_trt, site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species)%>%
  summarise(relcov=mean(relcov))%>%
  filter(exp_trt!='NIN::herbdiv::0::5F' & site_code!='GVN')
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
#   #merging back with labels to get back plot_mani
#   turnover=labels%>%
#     left_join(disappearance, by='calendar_year')%>%
#     left_join(appearance, by='calendar_year')%>%
#     filter(calendar_year==max(calendar_year))%>%
#     select(exp_trt, plot_mani, appearance, disappearance)
# 
#   #pasting variables into the dataframe made for this analysis
#   turnoverAll=rbind(turnover, turnoverAll)  
# }
# 
# turnoverCtl <- turnoverAll%>%
#   filter(plot_mani==0)%>%
#   separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)%>%
#   select(site_code, project_name, community_type, appearance, disappearance)
# names(turnoverCtl)[names(turnoverCtl)=='appearance'] <- 'appearance_ctl'
# names(turnoverCtl)[names(turnoverCtl)=='disappearance'] <- 'disappearance_ctl'
# 
# turnoverDiff <- turnoverAll%>%
#   mutate(trt=ifelse(plot_mani==0, 'ctl', 'trt'))%>%
#   separate(exp_trt, into=c('site_code', 'project_name', 'community_type', 'treatment'), sep='::', remove=F)%>%
#   filter(trt!='ctl')%>%
#   left_join(turnoverCtl, by=c('site_code', 'project_name', 'community_type'))%>%
#   mutate(appearance_diff=appearance-appearance_ctl, disappearance_diff=disappearance-disappearance_ctl)
# 
# # plot(turnoverDiff$plot_mani, turnoverDiff$appearance_diff)
# # plot(turnoverDiff$plot_mani, turnoverDiff$disappearance_diff)
# 
# turnoverRichness <- richness4%>%
#   left_join(turnoverDiff, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
#   select(site_code, project_name, community_type, treatment, experiment_length, plot_mani, intercept, slope, quad, min_year, nutrients, water, carbon, precip, alt_length, final_year_estimate, yr20, appearance_diff, disappearance_diff)%>%
#   filter(slope<0, quad>0)
# 
# plot(turnoverRichness$quad, turnoverRichness$appearance_diff)
# plot(turnoverRichness$quad, turnoverRichness$disappearance_diff)
# plot(turnoverRichness$final_year_estimate, turnoverRichness$appearance_diff)
# plot(turnoverRichness$final_year_estimate, turnoverRichness$disappearance_diff)
# plot(turnoverRichness$yr20, turnoverRichness$appearance_diff)
# plot(turnoverRichness$yr20, turnoverRichness$disappearance_diff)
# 
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
#   select(site_code, project_name, community_type, treatment, plot_mani, genus_species, first, last, experiment_length, intercept, slope, quad, nutrients, water, carbon, precip, alt_length, yr10, final_year_estimate)%>%
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

dispersionReps <- dispersion4%>%
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

ggplot(data=dispersionReps, aes(x=rep_num, y=yr10)) +
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

meanCompare <- mean4%>%
  left_join(expRawMean, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))


#plot without N at four factors, with N at five factors
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3), variable='final_year_estimate', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none')
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&herb_removal>0), variable='final_year_estimate', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none')
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plant_mani>0), variable='final_year_estimate', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
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











