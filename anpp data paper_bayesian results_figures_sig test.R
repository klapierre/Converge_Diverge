library(ggplot2)
library(grid)
library(mgcv)
library(plyr)
library(dplyr)
library(tidyr)

setwd('C:\\Users\\Kim\\Desktop\\bayesian output')
# setwd("~/Dropbox/Working GroupConverge Diverge/first data paper - bayesian results/anpp")


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

expInfo2 <- read.csv('ExperimentInformation_Mar2016.csv')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), precip=mean(precip))

rawData <- read.csv('ForBayesianAnalysisANPP_March2016b.csv')# get SD and means and use this backtransform the chains

mean(rawData$anpp_PC, na.rm=T) #0.3557605
sd(rawData$anpp_PC, na.rm=T) #0.5999784



################################################################################
################################################################################

###in our dropbox THIS IS FOR TREATMENT - Skinny lines for each treatment
chainsCommunity <- read.csv('fullChains_anpp.csv')

chainsMeanIntercept <- chainsCommunity[,8:182]%>%
  gather(key=parameter, value=value, B.1.1:B.175.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Intercept=median(value))%>%
  mutate(lower=Intercept-2*sd, upper=Intercept+2*sd, lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Intercepts=ifelse(diff==-2, 0, Intercept))%>%
  select(parameter, Intercepts)
names(chainsMeanIntercept)[1]<-"name"

chainsMeanSlope <- chainsCommunity[,183:357]%>%
  gather(key=parameter, value=value, B.1.2:B.175.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Slope=median(value))%>%
  mutate(lower=Slope-2*sd, upper=Slope+2*sd,lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Slopes=ifelse(diff==-2, 0, Slope))%>%
  select(parameter, Slopes)
names(chainsMeanSlope)[1]<-"name2"

chainsMeanQuad <- chainsCommunity[,358:532]%>%
  gather(key=parameter, value=value, B.1.3:B.175.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), Quad=median(value))%>%
  mutate(lower=Quad-2*sd, upper=Quad+2*sd,lower1=sign(lower), upper1=sign(upper), diff=lower1-upper1, Quads=ifelse(diff==-2, 0, Quad))%>%
  select(parameter, Quads)
names(chainsMeanQuad)[1]<-"name3"

#Select the interceipt, slope and quadratic columns

#Getting fat lines for each plot mani
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.1, U_int.2, U_int.3, U_int.4, U_lin.1, U_lin.2, U_lin.3, U_lin.4, U_quad.1, U_quad.2, U_quad.3, U_quad.4, mu_int, mu_lin, mu_quad)%>%
  gather(key=parameter, value=value, U_int.1:mu_quad)%>%
  group_by(parameter)%>%
  summarise(median=median(value))

###mean change
means<-cbind(chainsMeanQuad, chainsMeanSlope, chainsMeanIntercept)%>%
  select(-name2, -name3)%>%
  separate(name, into=c("head","exp","type"), sep="\\.", remove=F)%>%
  select(-head, -type)%>%
  mutate(exp2=as.integer(exp))

means2<-means[order(means$exp2),]

means3<-cbind(expInfo, means2)%>%
  select(-exp2, -exp)

mean<-merge(means3, expInfo2, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  #get estimates at 10 years
  mutate(final_year_estimate=(Intercepts + alt_length*Slopes + (alt_length^2)*Quads)*0.5999784 + 0.3557605)%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

#write.csv(mean, "to_clean.csv", row.names = F)

#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-2,4))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  ylim(-3,6) +
  xlab('Standardized Year') +
  ylab('ANPP Percent Change')

meanPlot <- meanPlot + 
#below are the individual treatment lines
stat_function(fun=function(x){(0.6226485 + 0*x + 0.02443635*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2915235*x + 0.0266958*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.493569 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.5123575 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.581095 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.608367 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.613646 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6061535 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6122845 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6090465 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7669475 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7489 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.821094 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.718188 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5120765 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(2.013385 + -0.1899725*x + 0.00625164*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.4830435 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(2.026155 + -0.177376*x + 0.01303065*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.6100635 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(1.7246 + -0.161294*x + 0.007621205*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(2.298985 + -0.157134*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.201294*x + 0.01700125*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(2.195625 + -0.462121*x + 0.026873*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.180092*x + 0.0152469*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(2.1797 + -0.4503375*x + 0.026079*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(1.985855 + -0.305241*x + 0.0186478*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5870305 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6883055 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.7468615 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6795425 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5306735 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.575866 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.555526 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7954255 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.8076045 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.835492 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.468859 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5418745 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.00554483*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.00625023*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.014045 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.28136 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2072205*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.343163*x + -0.011665*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.297059*x + -0.01111725*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.00657365*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1383565*x + -0.00818273*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.215171*x + -0.00950007*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.555906 + 0.1564105*x + -0.00857192*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(-0.646482 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.634151 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1330845*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4565855 + 0.106005*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6302455 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.701547 + 0.228308*x + -0.0176099*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6453665 + 0*x + -0.0147007*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7886345 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2183975*x + -0.01176735*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.762092 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.0288036*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.7820515 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.814463 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.650522 + 0*x + 0.0209688*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8356555 + 0*x + 0.02176765*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.02275205*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.06018 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9923395 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7667365 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.991824 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.01001 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.716893 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.831827 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.949299 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.10439 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.928593 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(1.036745 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.7487005 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7827025 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4643045 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4843895 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-0.743907 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6372555 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.514331 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5431495 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5335985 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=0.5, xlim=c(0,2), colour='#1400E544') +

  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.321176 + 0*x + -0*x^2)*0.5999784 + 0.3557605}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.321176 + 0.3400775) + 0*x + 0*x^2)*0.5999784 + 0.3557605}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.321176000 + 0*x + -0*x^2)*0.5999784 + 0.3557605}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.321176000 + 0*x + (0.001259430 +0.01104675)*x^2)*0.5999784 + 0.3557605}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.321176 + 2.30564) + (0.01842905 + -0.325715)*x + (0.00125943+0.0199754)*x^2)*0.5999784 + 0.3557605}, size=3, xlim=c(0,22), colour='#EC1804')

print(meanPlot)
#export at 1200x1000


###by resource mani

#mean change
meanResourceDrought <- mean%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
meanResource <- mean%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(meanResourceDrought)
  
# meanResourcePlot <- ggplot(data=barGraphStats(data=meanResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.5)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')

#dispersion change

dispersionResourceDrought <- dispersion%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
dispersionResource <- dispersion%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(dispersionResourceDrought)

# dispersionResourcePlot <- ggplot(data=barGraphStats(data=dispersionResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.08, 0.15)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.14, label='(b)', size=10, hjust='left')

#richness change

richnessResourceDrought <- richness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
richnessResource <- richness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(richnessResourceDrought)

# richnessResourcePlot <- ggplot(data=barGraphStats(data=richnessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.22, 0.16)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')

#evenness change
evennessResourceDrought <- evenness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated<0)%>%
  mutate(resource='drought')
evennessResource <- evenness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, carbon, precip, yr10, yr20, final_year_estimate)%>%
  gather(key=resource, value=manipulated, nutrients:precip)%>%
  filter(manipulated>0)%>%
  rbind(evennessResourceDrought)

# evennessResourcePlot <- ggplot(data=barGraphStats(data=evennessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.06)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left')
  
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))


#by resource at final year
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=meanResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(0, 0.5)) +
  xlab('')+
  annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')
dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=dispersionResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.08, 0.15)) +
  xlab('')+
  annotate('text', x=0.5, y=0.14, label='(b)', size=10, hjust='left')
richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=richnessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.22, 0.16)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')
evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=evennessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(0, 0.06)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))

# #20 year estimate
# meanResourcePlot20 <- ggplot(data=barGraphStats(data=meanResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.5)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')
# dispersionResourcePlot20 <- ggplot(data=barGraphStats(data=dispersionResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.08, 0.15, 0.04), name='Change in Dispersion') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.08, 0.15)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.14, label='(b)', size=10, hjust='left')
# richnessResourcePlot20 <- ggplot(data=barGraphStats(data=richnessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(-0.22, 0.16)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')
# evennessResourcePlot20 <- ggplot(data=barGraphStats(data=evennessResource, variable='final_year_estimate', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
#   scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
#                    labels=c('+CO2', '+nuts', '+rain', '-rain')) +
#   coord_cartesian(ylim=c(0, 0.06)) +
#   xlab('Resource Manipulated')+
#   annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlot20, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlot20, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlot20, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlot20, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))




#by resource mani boxplot
meanResourceBoxFinal <- ggplot(data=meanResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(0, 1.0)) +
  xlab('')+
  annotate('text', x=0.5, y=0.98, label='(a)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))
dispersionResourceBoxFinal <- ggplot(data=dispersionResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.1), name='Dispersion Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.4, 0.4)) +
  xlab('')+
  annotate('text', x=0.5, y=0.39, label='(b)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))
richnessResourceBoxFinal <- ggplot(data=richnessResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-1, 0.6, 0.3), name='Proportion Richness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-1, 0.6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.58, label='(c)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))
evennessResourceBoxFinal <- ggplot(data=evennessResource, aes(x=resource, y=final_year_estimate)) +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-0.4, 0.6, 0.2), name='Evenness Change') +
  scale_x_discrete(limits=c('carbon', 'nutrients', 'drought', 'precip'),
                   labels=c('+CO2', '+nuts', '+rain', '-rain')) +
  coord_cartesian(ylim=c(-0.4, 0.6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.58, label='(d)', size=10, hjust='left') + 
  geom_hline(aes(yintercept=0))


pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourceBoxFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourceBoxFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourceBoxFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourceBoxFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))



#other inset options

#inset - density plot of mean change (all datapoints)
mean10yr <- ggplot(data=mean, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
mean20yr <- ggplot(data=mean, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
meanFinalYear <- ggplot(data=mean, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')
dispersion10yr <- ggplot(data=dispersion, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
dispersion20yr <- ggplot(data=dispersion, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
dispersionFinalYear <- ggplot(data=dispersion, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')
richness10yr <- ggplot(data=richness, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
richness20yr <- ggplot(data=richness, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
richnessFinalYear <- ggplot(data=richness, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')
evenness10yr <- ggplot(data=evenness, aes(x=yr10)) +
  geom_density() +
  xlab('') +
  ylab('')
evenness20yr <- ggplot(data=evenness, aes(x=yr20)) +
  geom_density() +
  xlab('') +
  ylab('')
evennessFinalYear <- ggplot(data=evenness, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('') +
  ylab('')


pushViewport(viewport(layout=grid.layout(3,4)))
print(mean10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(mean20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(meanFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 1))
print(dispersion10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(dispersion20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(dispersionFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 2))
print(richness10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(richness20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
print(richnessFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 3))
print(evenness10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(evenness20yr, vp=viewport(layout.pos.row = 2, layout.pos.col = 4))
print(evennessFinalYear, vp=viewport(layout.pos.row = 3, layout.pos.col = 4))

      
      
      
###summary stats from bayesian output
chainsCommunitySummary <- chainsCommunity%>%
  select(U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4, E_int.1.1, E_int.2.1, E_int.3.1, E_int.4.1, E_slope.1.1, E_slope.2.1, E_slope.3.1, E_slope.4.1, E_quad.1.1, E_quad.2.1, E_quad.3.1, E_quad.4.1, E_int.1.2, E_int.2.2, E_int.3.2, E_int.4.2, E_slope.1.2, E_slope.2.2, E_slope.3.2, E_slope.4.2, E_quad.1.2, E_quad.2.2, E_quad.3.2, E_quad.4.2, D_int.1.1, D_int.2.1, D_int.3.1, D_int.4.1, D_slope.1.1, D_slope.2.1, D_slope.3.1, D_slope.4.1, D_quad.1.1, D_quad.2.1, D_quad.3.1, D_quad.4.1, D_int.1.2, D_int.2.2, D_int.3.2, D_int.4.2, D_slope.1.2, D_slope.2.2, D_slope.3.2, D_slope.4.2, D_quad.1.2, D_quad.2.2, D_quad.3.2, D_quad.4.2)%>%
  gather(key=parameter, value=value, U_int.2.1:D_quad.4.2)%>%
  group_by(parameter)%>%
  summarise(median=median(value), sd=sd(value))%>%
  mutate(CI=sd*2)

write.csv(chainsCommunitySummary, 'bayesian_summary stats.csv')
chainsCommunitySummary2 <- read.csv('bayesian_summary stats.csv')

#mean plots
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='mean change'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

#dispersion plots
dispersionIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

dispersionSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

dispersionQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

#richness plots
richnessIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

richnessSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

richnessQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.035, 0.05) +
  coord_flip()

#evenness plots
evennessIntPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='int'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

evennessSlopePlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='slope'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-0.5, 0.6) +
  coord_flip()

evennessQuadPlot <- ggplot(data=subset(chainsCommunitySummary2, parameter=='quad'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'gamma', 'mani 5', 'mani 4', 'mani 3', 'mani 2'),
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
meanOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='mean change'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='dispersion'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='richness'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary2, variable=='evenness'&predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('quad', 'slope', 'int'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

pushViewport(viewport(layout=grid.layout(1,4)))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
#export at 1800x1200
