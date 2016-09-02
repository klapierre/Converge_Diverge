library(ggplot2)
library(grid)
library(mgcv)
library(plyr)
library(dplyr)
library(tidyr)

setwd("C:\\Users\\Kim\\Desktop\\bayesian output\\9 yr subset\\anpp")


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
expInfo <- read.csv('ExperimentInformation_Mar2016_anpp.csv')

expInfo2 <- read.csv('C:\\Users\\Kim\\Desktop\\bayesian output\\ExperimentInformation_Mar2016.csv')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), precip=mean(precip))

rawData <- read.csv('ForBayesianAnalysisANPP_9yr_July2016.csv')# get SD and means and use this backtransform the chains

mean(rawData$anpp_PC, na.rm=T) #0.3275504
sd(rawData$anpp_PC, na.rm=T) #0.5319707



################################################################################
###############################################################################
#get model output
chains1 <- read.csv('results_anpp_9yr_0.csv', comment.char='#')
chains1 <- chains1[-1:-5000,]
chains2 <- read.csv('results_anpp_9yr_1.csv', comment.char='#')
chains2 <- chains2[-1:-5000,]
chains3 <- read.csv('results_anpp_9yr_2.csv', comment.char='#')
chains3 <- chains3[-1:-5000,]
chains4 <- read.csv('results_anpp_9yr_3.csv', comment.char='#')
chains4 <- chains4[-1:-5000,]

chainsCommunity <- rbind(chains1, chains2, chains3, chains4)

chainsMeanIntercept <- chainsCommunity[,7:177]%>%
  gather(key=parameter, value=value, B.1.1:B.171.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Intercept=median(value))%>%
  mutate(lower=Intercept-2*sd, upper=Intercept+2*sd, lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Intercepts=ifelse(diff==-2, 0, Intercept))%>%
  select(parameter, Intercepts)
names(chainsMeanIntercept)[1]<-"parameter1"

chainsMeanSlope <- chainsCommunity[,178:348]%>%
  gather(key=parameter, value=value, B.1.2:B.171.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Slope=median(value))%>%
  mutate(lower=Slope-2*sd, upper=Slope+2*sd,lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Slopes=ifelse(diff==-2, 0, Slope))%>%
  select(parameter, Slopes)
names(chainsMeanSlope)[1]<-"parameter2"

chainsMeanQuad <- chainsCommunity[,349:519]%>%
  gather(key=parameter, value=value, B.1.3:B.171.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Quad=median(value))%>%
  mutate(lower=Quad-2*sd, upper=Quad+2*sd,lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Quads=ifelse(diff==-2, 0, Quad))%>%
  select(parameter, Quads)
names(chainsMeanQuad)[1]<-"parameter3"

#Select the intercept, slope and quadratic columns

#Getting fat lines for each plot mani
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.1, U_int.2, U_int.3, U_int.4, U_lin.1, U_lin.2, U_lin.3, U_lin.4, U_quad.1, U_quad.2, U_quad.3, U_quad.4, mu_int, mu_lin, mu_quad)%>%
  gather(key=parameter, value=value, U_int.1:mu_quad)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))

###mean change
means<-cbind(chainsMeanQuad, chainsMeanSlope, chainsMeanIntercept)%>%
  select(-parameter2, -parameter3)%>%
  separate(parameter1, into=c("head","exp","type"), sep="\\.", remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

means2<-means[order(means$exp2),]
means3<-cbind(expInfo, means2)%>%
  select(-exp2, -exp)

mean<-left_join(means3, expInfo2, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  #get estimates at 10 years
  mutate(final_year_estimate=(Intercepts + alt_length*Slopes + (alt_length^2)*Quads)*0.5319707 + 0.3275504)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below


#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(xlim=c(0,9), ylim=c(-2,4))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # ylim(-3,6) +
  xlab('Standardized Year') +
  ylab('ANPP Percent Change')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(1.046935 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8064525 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.755492 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.7803395 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8964485 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.745116 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.853198 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7770785 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.744075 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.201915 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(2.13245 + 0*x + 0.05324905*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.601425 + 0*x + 0.057266*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.60898 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(3.76809 + -0.963509*x + 0.0919276*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(3.50692 + -1.109575*x + 0.0874495*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.4178585*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(4.03322 + -1.467985*x + 0.121428*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(2.75682 + -0.869359*x + 0.0887705*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.02261 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.016115 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0886347*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7162725 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8301865 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.9063555 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.821184 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9757185 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.244555 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.288215 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.519347*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(1.101905 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7908585 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.740348 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.4184 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00939 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.24107 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.182135 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.105705 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.45114 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.540485 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.363265 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.41712 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.884535 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.60223 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.977796 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.943845 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.9664935 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0591111*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0634869*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4127685*x + -0.07210295*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0666456*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.541249*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4764625*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(-0.8434805 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  

  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.32135300 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.32135300+0.35868200) + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.32135300 + (0.090636850-0.40954050)*x + (-0.014275050+0.08573735)*x^2)*0.5319707 + 0.3275504}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.32135300 + 0*x + 0*x^2)*0.5319707 + 0.3275504}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.32135300+2.84595000) + (0.090636850-0.78978150)*x + (-0.014275050+0.08199180)*x^2)*0.5319707 + 0.3275504}, size=3, xlim=c(0,8), colour='#EC1804')

print(meanPlot)
#export at 1200x1000





##############################################################################################
#absolute value of ANPP percent change (this is to see magnitude, regardless of direction, of response)

#raw data
rawDataa <- read.csv('ForBayesianAnalysisANPP_9 yr_abs value_ANPP_Aug2016.csv')# get SD and means and use this backtransform the chains

mean(rawDataa$anpp_PC, na.rm=T) #0.4028107
sd(rawDataa$anpp_PC, na.rm=T) #0.4512122

#get model output
chains1a <- read.csv('results_anpp_absval_9yr_0.csv', comment.char='#')
chains1a <- chains1a[-1:-5000,]
chains2a <- read.csv('results_anpp_absval_9yr_1.csv', comment.char='#')
chains2a <- chains2a[-1:-5000,]
chains3a <- read.csv('results_anpp_absval_9yr_2.csv', comment.char='#')
chains3a <- chains3a[-1:-5000,]
chains4a <- read.csv('results_anpp_absval_9yr_3.csv', comment.char='#')
chains4a <- chains4a[-1:-5000,]

chainsCommunitya <- rbind(chains1a, chains2a, chains3a, chains4a)

chainsMeanIntercepta <- chainsCommunitya[,7:180]%>%
  gather(key=parameter, value=value, B.1.1:B.174.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Intercept=median(value))%>%
  mutate(lower=Intercept-2*sd, upper=Intercept+2*sd, lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Intercepts=ifelse(diff==-2, 0, Intercept))%>%
  select(parameter, Intercepts)
names(chainsMeanIntercepta)[1]<-"parameter1"

chainsMeanSlopea <- chainsCommunitya[,181:354]%>%
  gather(key=parameter, value=value, B.1.2:B.174.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Slope=median(value))%>%
  mutate(lower=Slope-2*sd, upper=Slope+2*sd,lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Slopes=ifelse(diff==-2, 0, Slope))%>%
  select(parameter, Slopes)
names(chainsMeanSlopea)[1]<-"parameter2"

chainsMeanQuada <- chainsCommunitya[,355:528]%>%
  gather(key=parameter, value=value, B.1.3:B.174.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Quad=median(value))%>%
  mutate(lower=Quad-2*sd, upper=Quad+2*sd,lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Quads=ifelse(diff==-2, 0, Quad))%>%
  select(parameter, Quads)
names(chainsMeanQuada)[1]<-"parameter3"

#Select the intercept, slope and quadratic columns

#Getting fat lines for each plot mani
chainsCommunity2a <- chainsCommunitya%>%
  select(lp__, U_int.1, U_int.2, U_int.3, U_int.4, U_lin.1, U_lin.2, U_lin.3, U_lin.4, U_quad.1, U_quad.2, U_quad.3, U_quad.4, mu_int, mu_lin, mu_quad)%>%
  gather(key=parameter, value=value, U_int.1:mu_quad)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))

###mean change
meansa<-cbind(chainsMeanQuada, chainsMeanSlopea, chainsMeanIntercepta)%>%
  select(-parameter2, -parameter3)%>%
  separate(parameter1, into=c("head","exp","type"), sep="\\.", remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

means2a<-meansa[order(meansa$exp2),]
means3a<-cbind(expInfo, means2a)%>%
  select(-exp2, -exp)

meana<-left_join(means3a, expInfo2, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  #get estimates at 10 years
  mutate(final_year_estimate=(Intercepts + alt_length*Slopes + (alt_length^2)*Quads)*0.4512122 + 0.4028107)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below


#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(xlim=c(0,9), ylim=c(-2,4))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # ylim(-3,6) +
  xlab('Standardized Year') +
  ylab('ANPP Percent Change')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(1.05828 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.8281135 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8430505 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8367905 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.991647 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.773689 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7087565 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.757687 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.719149 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.87866 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(1.262895 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(2.29336 + -0.4342045*x + 0.05992285*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.63343 + 0*x + 0.0589198*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(1.69643 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(4.30099 + -1.138175*x + 0.1039075*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(3.89596 + -1.169645*x + 0.0957387*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(4.55578 + -1.70178*x + 0.1396845*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(3.090745 + -1.032235*x + 0.105121*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.9976815 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.952915 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0676505*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.722763 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(1.064865 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0858963*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0915937*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0.901087 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.514397*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0928288*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(1.45057 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7249785 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7688775 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.522345 + 0*x + 0.0923909*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.84412 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(2.07543 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.9062275 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7448915 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0618842*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.05704595*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.53897*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
  
  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.27625100 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=3, xlim=c(0,8), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.27625100+0.41339050) + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=3, xlim=c(0,8), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.27625100 + 0*x + (-0.009745275+0.08072975)*x^2)*0.4512122 + 0.4028107}, size=3, xlim=c(0,8), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.27625100 + 0*x + 0*x^2)*0.4512122 + 0.4028107}, size=3, xlim=c(0,8), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.27625100+3.09588000) + (0.058578050-0.87264700)*x + (-0.009745275+0.08948375)*x^2)*0.4512122 + 0.4028107}, size=3, xlim=c(0,8), colour='#EC1804')

print(meanPlot)
#export at 1200x1000

