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

rawData <- read.csv('9 yr subset\\ForBayesianAnalysis_9yr_July2016.csv')

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
chains1 <- read.csv('9 yr subset\\diversity_MV_rdisp_9yr_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('9 yr subset\\diversity_MV_rdisp_9yr_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('9 yr subset\\diversity_MV_rdisp_9yr_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('9 yr subset\\diversity_MV_rdisp_9yr_3.csv', comment.char='#')
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
chainsMeanIntercept <- chainsCommunity[,7:292]%>%
  gather(key=parameter, value=value, B.1.1.1:B.286.1.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsMeanIntercept)[1] <- 'parameter1'

chainsMeanSlope <- chainsCommunity[,1151:1436]%>%
  gather(key=parameter, value=value, B.1.1.2:B.286.1.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsMeanSlope)[1] <- 'parameter2'

chainsMeanQuad <- chainsCommunity[,2295:2580]%>%
  gather(key=parameter, value=value, B.1.1.3:B.286.1.3)%>%
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
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.1467006+0.2945762,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.1467006+0.2945762,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.1467006+0.2945762)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,',
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
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.624931 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.497628 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4608355 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.262779*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.268465*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.6839465 + 0.304835*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.417785 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.500515 + 0.310005*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.41333 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.3392 + 0.3550595*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.615309 + 0.2829645*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.76686 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.6198035 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.476617 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.773173 + 0.364923*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5139865 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.280995*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.324964*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.896618 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7948745 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.136515 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.842226 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-0.615325 + 0.306639*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.898885 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9184725 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-1.705895 + 0.320535*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.8787 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.830485 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.902805 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.753445 + 0.3226425*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.867645 + 0.288233*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5735335 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4719225 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.570704 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6043165 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4978665 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.325028*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6263375 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3506*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6515585 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.819825 + 0.3696245*x + -0.04884075*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8221565 + 0.381263*x + -0.0426214*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(-1.35235 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0.5310735 + 0.637189*x + -0.0503678*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-1.07008 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.500106 + 0.6688475*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.601428 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.7568625*x + -0.0546293*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(1.10891 + 0.7460085*x + -0.04506775*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5886435 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4778975*x + -0.03700415*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(1.015045 + 0.8836565*x + -0.05745245*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5972225 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.731278*x + -0.04109165*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7967945 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.60624 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8183435 + 0.3585765*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.137325 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-1.443055 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.32975 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3291275*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.927291 + 0.302627*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3292065*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7933105 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.557132 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.785712 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.845585 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8883225 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.216325 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.05126 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.105115 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.14704 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9870105 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5791375 + 0.25823*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4383955 + 0.2390025*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.497657 + 0.327566*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.470957 + 0.261182*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.5032405 + 0.327268*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5194505 + 0.2882635*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.373941*x + -0.03088135*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5044105 + 0.239516*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5140035 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4269215 + 0.2968065*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6334475 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5986235 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6884565 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.551025 + 0.2791055*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.526847 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(2.100915 + 0.351319*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8365345 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.698142 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.921241 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3892825*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3562725*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3484035*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.829622 + 0.306261*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7713635 + 0.383259*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(-0.754456 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.798626 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6594935 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.840213 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.055115 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.14343 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00193 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.916156 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8668245 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6920445 + 0.403295*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.506316 + 0.442664*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9017265 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.089005 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.932829 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.057795 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6079395 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.251449*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.5088685 + 0.2684395*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4123305*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.435804*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.2735075*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.305237*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.821825 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7867905 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.719899 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.8358705 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.476713 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.465977 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.455085 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.874279 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.593283 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.889418 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.556552 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.00579 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6194595 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9423925 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5805745 + 0.2932945*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.598302 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.54133 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.329631*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.489523 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.89324 + 0.389329*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.377559*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.567809 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.757067 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4027395*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3866605*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.357521*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.424607*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(-0.834602 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.23529 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.729561 + 0.3994485*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.21894 + 0.259638*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.44082 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.180855 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.246335 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.097335 + 0.253746*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-1.08749 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.933388 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9801255 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.905735 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8215465 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.38898 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7941575 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6121335 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.624009 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.780843 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.841509 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5494765 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7121945 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7142395 + 0.380947*x + -0.05612305*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.484405 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.40595 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.08429 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.07238 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.1046 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.0753 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.652563 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8329265 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7584335 + 0.306344*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8889505 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.519804 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.303634*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.937693 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.548771 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.958595 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.711629 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6805305 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.355481*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.4050825*x + -0.05715535*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9750215 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.571081 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0.7702325 + 0*x + 0*x^2)*(0.1467006)+(0.2945762)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.2945762)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0.2945762)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.51288850 + 0.18689450*x + 0*x^2)*0.1467006 + 0.2945762}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.51288850+0.16444100) + 0.18689450*x + 0*x^2)*0.1467006 + 0.2945762}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.51288850 + 0.18689450*x + 0*x^2)*0.1467006 + 0.2945762}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.51288850 + 0.18689450*x + 0*x^2)*0.1467006 + 0.2945762}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.51288850+0.62378000) + (0.18689450+0.41136250)*x + 0*x^2)*0.1467006 + 0.2945762}, size=3, xlim=c(0,22), colour='#EC1804')

# print(meanPlot) #export at 1200x1000




###dispersion
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsDispersionIntercept <- chainsCommunity[,293:578]%>%
  gather(key=parameter, value=value, B.1.2.1:B.286.2.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsDispersionIntercept)[1] <- 'parameter1'

chainsDispersionSlope <- chainsCommunity[,1437:1722]%>%
  gather(key=parameter, value=value, B.1.2.2:B.286.2.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsDispersionSlope)[1] <- 'parameter2'

chainsDispersionQuad <- chainsCommunity[,2581:2866]%>%
  gather(key=parameter, value=value, B.1.2.3:B.286.2.3)%>%
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
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.08810367-0.000327455,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.08810367-0.000327455,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.08810367-0.000327455)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, 'grey', ifelse(plot_mani==2, 'grey', ifelse(plot_mani==3, 'grey', ifelse(plot_mani==4, 'grey', 'grey')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(dispersion4,'dispersion_equations.csv', row.names=F)



#main figure
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.5))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.5, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.646359 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.7279545 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.05558085*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.08668 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5101515*x + 0.07495585*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,15), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.866 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.05157295*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7529325 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0558274*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0529491*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.762412 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(-2.01068 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4450005*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0.503992*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0.354346*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,20), colour='grey') +
  stat_function(fun=function(x){(0 + 0.332514*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,20), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.63161 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.333575*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4242495*x + -0.04899575*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.774274 + -0.54049*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.12988 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.7966605 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6446805 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03273815*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03239575*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04045715*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03714*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.03282775*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0.26028*x + -0.0400182*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(-2.41455 + 0.8363315*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.33911 + 0.723073*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.80573 + 0.604775*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.40907 + 0.620335*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.745865 + 0.747064*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.32393 + 0.679369*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.24024 + 0.865218*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.64644 + 0.817315*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7183055 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.049055 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.726168 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.7145155 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.848956 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.834945 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.80645 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(-0.794913 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.278835 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04864775*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.936478 + 0.598514*x + -0.07231415*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.801448 + 0.52668*x + -0.07019735*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.04556425*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0467771*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4181375*x + -0.0797539*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.931164 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.7559425 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.501099*x + 0.06824365*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(1.028255 + -0.5410785*x + 0.05955305*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.538687*x + 0.0627168*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.07898115*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4340205*x + 0.0905267*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.5505165*x + 0.1009525*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4714485*x + 0.0731279*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3537435*x + 0.05269045*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.397942*x + 0.04994215*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0377521*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.344031*x + 0.05539215*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.343865*x + 0.05269605*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9584225 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.741638 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.941892 + 0*x + 0.0447516*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.01752 + 0*x + 0.0487451*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.25957 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.7258245 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.13249 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.664125 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.55072 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7910595 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.671217 + 0.4156465*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.9686245 + 0.4474155*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0503415*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.647347 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(-0.6578385 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.760824 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9227695 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.8629865 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.14327 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.4297565*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.004965 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.3735735*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.0567558*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.05685795*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.707714 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.681632 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.147725 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.6390085 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.8152845 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.822044 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.393605 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.9126305 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.67917 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.483345 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.886816 + 0*x + 0*x^2)*(0.08810367)+(-0.000327455)}, size=0.5, xlim=c(0,2), colour='grey') +
#estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0 - 0}, size=3, xlim=c(0,23), colour='black')


# print(dispersionPlot) #export at 1200x1000







###richness
chainsRichnessIntercept <- chainsCommunity[,865:1150]%>%
  gather(key=parameter, value=value, B.1.4.1:B.286.4.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsRichnessIntercept)[1] <- 'parameter1'

chainsRichnessSlope <- chainsCommunity[,2009:2294]%>%
  gather(key=parameter, value=value, B.1.4.2:B.286.4.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsRichnessSlope)[1] <- 'parameter2'

chainsRichnessQuad <- chainsCommunity[,3153:3438]%>%
  gather(key=parameter, value=value, B.1.4.3:B.286.4.3)%>%
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
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.2160529-0.05773279,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.2160529-0.05773279,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.2160529-0.05773279)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(richness4,'richness_equations.csv', row.names=F)



#main figure
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1,0.8))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6678115 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.698509 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.056561*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.367244*x + 0.0616186*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.0503544*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0599163*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4288385*x + 0.06648325*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5613775 + -0.431911*x + 0.0688692*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5399525 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3523535*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5490425 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.923374 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.651411 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8512955 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.950261 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.659194 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.399205*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6834495 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0.5902645 + -0.398915*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.1023*x + 0.0915164*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.9938925 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.913884*x + 0.0726331*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.9307525*x + 0.081473*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.803972 + -0.6712745*x + 0.0787298*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.3443705*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.157755*x + 0.102515*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -1.13361*x + 0.0972988*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.942401*x + 0.089146*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0.7596585 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.343147*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.17964 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.6595975 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.5364055*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.6580555 + -0.46798*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.580661 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.3927705*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0.9011705 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.9396295 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.549832 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.03477885*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3277205*x + -0.0402905*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.266333*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.260388*x + -0.0331967*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.308783*x + -0.0370216*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5921705 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.544612 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.5960885 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4049965*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4455545*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5579015 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7988 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.682074 + 0*x + -0.050763*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.73344 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.05297675*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.999837 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.04817 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.9214395 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.06124 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.649476 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(1.004075 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.5616845 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6090475 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.3031305*x + 0.05366395*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.04379795*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.73553 + -0.397937*x + 0.0526999*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.4075325*x + 0.05244945*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.04983795*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.322251*x + 0.04497575*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.04687745*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.5811615 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.154385 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.565645 + 0.5295125*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.9716425 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.505245 + 0.6704645*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4006115*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.478293*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.8453725*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6471325 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.905427 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.981171 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.720884 + 0*x + -0.0469464*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0.724636 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6902445 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.5635675 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.5794115 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.667049 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0559458*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4382435*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.8440905 + 0.383769*x + -0.0479365*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0.6544465 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.590637 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7390145 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.687125 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7353465 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.62558 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.849787 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.594627 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.649118 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.851535 + -0.4815675*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(1.180065 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.7091235 + -0.4389815*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.717867 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4498955*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.588809 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-1.156715 + 0*x + 0*x^2)*(0.2160529)+(-0.05773279)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(0.23707350 + -0.09856190*x + 0*x^2)*0.2160529 - 0.05773279}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(0.23707350 + -0.09856190*x + -0.01907205*x^2)*0.2160529 - 0.05773279}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(0.23707350 + -0.09856190*x + 0*x^2)*0.2160529 - 0.05773279}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(0.23707350 + -0.09856190*x + 0*x^2)*0.2160529 - 0.05773279}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(0.23707350 + (-0.09856190-0.77843500)*x + (0.0126042580+0.06362400)*x^2)*0.2160529 - 0.05773279}, size=3, xlim=c(0,22), colour='#EC1804')

# print(richnessPlot) #export at 1200x1000






###evenness
chainsEvennessIntercept <- chainsCommunity[,579:864]%>%
  gather(key=parameter, value=value, B.1.3.1:B.286.3.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-2*sd, upper=intercept+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsEvennessIntercept)[1] <- 'parameter1'

chainsEvennessSlope <- chainsCommunity[,1723:2008]%>%
  gather(key=parameter, value=value, B.1.3.2:B.286.3.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-2*sd, upper=slope+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsEvennessSlope)[1] <- 'parameter2'

chainsEvennessQuad <- chainsCommunity[,2867:3152]%>%
  gather(key=parameter, value=value, B.1.3.3:B.286.3.3)%>%
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
  mutate(yr9=(intercept + 8*slope + (10^2)*quad)*0.1005463+0.01731845,
         yr20=(intercept + 20*slope + (20^2)*quad)*0.1005463+0.01731845,
         final_year_estimate=(intercept + alt_length*slope + (alt_length^2)*quad)*0.1005463+0.01731845)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, intercept, curve2, slope, curve3, quad, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
write.csv(evenness4,'evenness_equations.csv', row.names=F)


#main figure
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.35,0.6))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.5253345*x + -0.05141035*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.674692 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6188155 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.418989*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.05712325*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.04888465*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6167815 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0365131*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.715968 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7809065 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.006755 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.789567 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.891855 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.04923835*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0510089*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4113055*x + -0.056961*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.049717*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.403136*x + -0.0428713*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.475638*x + -0.05203435*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.484005*x + -0.114774*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.7611505*x + -0.0529416*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.5681935*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.631419*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 1.75109*x + -0.150164*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.803356 + 1.063595*x + -0.07795585*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.78682*x + -0.06287555*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.594922 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.603927 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.784603 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.763337 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.747204 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.636035 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.617139 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.55719 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.611048 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.36634*x + -0.0387433*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.348132*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.68313 + 0.5000675*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.57599 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.41703 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.707962*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.800938*x + -0.06373195*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 1.33646*x + -0.0987586*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.401526*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.486728*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4378115*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4583705*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.426876*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.5491825*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.603334*x + -0.04966595*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.3449805*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.444089*x + -0.0434693*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.722594*x + -0.05868105*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.391688*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4571105*x + 0.0495247*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.623102 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.61407 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.636807 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5875515 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6008855 + 0.3949645*x + -0.04725425*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5842405 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.711308 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.407567*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.537223*x + 0.0645293*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4026695*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6295385 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9990005 + 0.839517*x + -0.06700885*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(-1.01303 + 0.7592195*x + -0.0603822*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.428663*x + -0.0558694*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.04956425*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1005463)+(0.01731845)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.21478250 + 0.10652350*x + -0.01704330*x^2)*0.1005463 + 0.01731845}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.21478250 + 0.10652350*x + -0.01704330*x^2)*0.1005463 + 0.01731845}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.21478250 + 0.10652350*x + -0.01704330*x^2)*0.1005463 + 0.01731845}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.21478250 + 0.10652350*x + -0.01704330*x^2)*0.1005463 + 0.01731845}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){(-0.21478250 + (0.10652350+0.84695100)*x + (-0.01704330-0.05836085)*x^2)*0.1005463 + 0.01731845}, size=3, xlim=c(0,22), colour='#EC1804')

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
  coord_cartesian(ylim=c(-0.08, 0.08)) +
  xlab('')+
  annotate('text', x=0.5, y=0.038, label='(b)', size=10, hjust='left')
richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=richnessResource, variable='yr9', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+' ~CO[2], '+nutrients', '+' ~H[2]*O, '-' ~H[2]*O)) +
  coord_cartesian(ylim=c(-0.11, 0.1)) +
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
  scale_y_continuous(limits=c(-0.2, 0.5), breaks=seq(-0.4, 0.4, 0.4)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Proportion\nRichness Change') +
  annotate('text', x=3.45, y=-0.2, label='(c)', size=10, hjust='left')

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

meanCompare <- mean4%>%
  left_join(expRawMean, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))


#plot without N at four factors, with N at five factors
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none')
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&herb_removal>0), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none')
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plant_mani>0), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
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











