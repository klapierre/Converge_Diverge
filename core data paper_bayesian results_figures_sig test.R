library(ggplot2)
library(grid)
library(mgcv)
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
expInfo <- read.csv('ExperimentInformation_Mar2016.csv')%>%
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

#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsMeanIntercept <- chainsCommunity[,7:296]%>%
  gather(key=parameter, value=value, B.1.1.1:B.290.1.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-sd, upper=intercept+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsMeanIntercept)[1] <- 'parameter1'

chainsMeanSlope <- chainsCommunity[,1167:1456]%>%
  gather(key=parameter, value=value, B.1.1.2:B.290.1.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-sd, upper=slope+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsMeanSlope)[1] <- 'parameter2'

chainsMeanQuad <- chainsCommunity[,2327:2616]%>%
  gather(key=parameter, value=value, B.1.1.3:B.290.1.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-sd, upper=quad+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
  select(parameter, quad)
names(chainsMeanQuad)[1] <- 'parameter3'

#get values for overall (mean) lines across levels of plot mani
#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4)%>%
  gather(key=parameter, value=value, U_int.2.1:mu_quad.4)%>%
  group_by(parameter)%>%
  summarise(median=median(value))

###mean change
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
  mutate(yr10=intercept + 10*slope + (10^2)*quad,
         yr20=intercept + 20*slope + (20^2)*quad,
         final_year_estimate=intercept + alt_length*slope + (alt_length^2)*quad)%>%
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

# #change in mean at 10 years vs. gamma diversity
# ggplot(data=mean, aes(x=rrich, y=yr10)) +
#   geom_point(aes(color=as.factor(plot_mani))) +
#   stat_smooth(method='gam', se=T, formula = y ~ s(x), color='black') +
#   scale_color_manual(values=c('#1400E5', '#4A06AC', '#800C74', '#B6123C', '#EC1804')) +
#   xlab('Gamma Diversity') +
#   ylab('Mean Change at 10 Years') +
#   scale_y_continuous(limits=c(0,1)) +
#   scale_x_continuous(breaks=seq(0,140,20)) +
#   theme(legend.position='none')


# #inset - density plot of mean change (all datapoints)
# meanInset <- ggplot(data=mean, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Mean Change') +
#   ylab('Density')


#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  ylim(-1,2) +
  xlab('Standardized Year') +
  ylab('Mean Change')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){(-0.4883495 + 0.126389*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.304715 + 0.141138*x + -0.00896859*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.292111 + 0.180991*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.631301 + 0.1862475*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5098295 + 0.1263765*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.2581595 + 0.1643195*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.3851125 + 0.1716785*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.331873 + 0.17049*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2023475*x + -0.00695669*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3590435 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.314034 + 0.2527795*x + -0.00880028*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5530665 + 0.3011225*x + -0.00937814*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.3663745 + 0.331224*x + -0.01298095*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.7045975 + 0.383736*x + -0.015702*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.27143 + 0.318652*x + -0.0198106*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.12772 + 0*x + -0.01042595*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.02848 + 0*x + -0.01078035*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.00456 + 0.117365*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.42911 + 0.2339005*x + -0.0123622*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.577968 + 0.210579*x + -0.01165585*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.4461505 + 0.206251*x + -0.01132075*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.333708 + 0.1469165*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.248163*x + -0.0116109*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.582217 + 0.266819*x + -0.0141823*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.3638485 + 0.2128925*x + -0.0115557*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1917215*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1748035*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1837935*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.18206*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2176955 + 0.2203325*x + -0.00871936*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8106405 + 0.1352605*x + -0.008917425*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.683987 + 0.09437635*x + -0.007811345*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.02411 + 0.1091545*x + -0.007566935*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7965215 + 0.202106*x + -0.0113552*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-0.609518 + 0.260602*x + -0.0151973*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8540235 + 0.1386425*x + -0.00978064*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8374905 + 0.07447665*x + -0.0062092*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-1.51362 + 0.184485*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.63217 + 0.126566*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.59007 + 0.0969311*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.63067 + 0.1041025*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.50476 + 0.148761*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-1.662305 + 0.1694735*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6354955 + 0*x + -0.01201635*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5043855 + 0*x + -0.01162175*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1745365*x + -0.01733725*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1469635*x + -0.0160752*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1145625*x + -0.0146415*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5628235 + 0*x + -0.01309805*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4731355 + 0.1710675*x + -0.017348*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6129605 + 0.09382325*x + -0.0142133*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4940195 + 0.0904691*x + -0.0144598*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.3754445 + 0*x + -0.01330685*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.0993544*x + -0.01370555*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.340781 + 0*x + -0.0112854*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1522415*x + -0.0156804*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1837785*x + -0.01724645*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.429411 + 0.08408785*x + -0.0126735*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.396701 + 0.1867935*x + -0.01425125*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6138175 + 0*x + -0.009181625*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.2374215 + 0.1762135*x + -0.0142518*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2457725*x + -0.0170403*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.669725 + 0.103378*x + -0.0107635*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.671122 + 0.08961935*x + -0.00680747*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.707299 + 0.1289235*x + -0.00888552*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(-1.30051 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
                stat_function(fun=function(x){(0.589338 + 0.3214405*x + -0.01045525*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.8919455 + 0.0441775*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.445834 + 0.530661*x + -0.01611715*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.4873165 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
                stat_function(fun=function(x){(0.6912785 + 0.2038135*x + -0.00703984*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.1335225*x + -0.008427285*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
                stat_function(fun=function(x){(1.248975 + 0.3036995*x + -0.01604555*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5986995 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
                stat_function(fun=function(x){(0 + 0.2480395*x + -0.0094451*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.229657 + 0*x + -0.006628955*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
                stat_function(fun=function(x){(1.026485 + 0.543711*x + -0.0258208*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.5274785 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
                stat_function(fun=function(x){(0.2401285 + 0.4682845*x + -0.0170641*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.7563585 + 0.109698*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.57131 + 0.139115*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7582515 + 0.2581355*x + -0.01034355*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.05297 + 0.172808*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-1.348815 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.24015 + 0.1212755*x + -0.0104831*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.19757*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4165555 + 0.176853*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4400815 + 0.3280515*x + -0.00898988*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.844727 + 0.2267175*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.484147 + 0.269437*x + -0.00988593*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
                stat_function(fun=function(x){(0 + 0.222553*x + -0.008679325*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.197577*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7693165 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.544599 + 0.182794*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(-0.63166 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.80487 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.85503 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.902496 + 0.129825*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-1.193055 + 0.113225*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.014385 + 0.1534395*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.103805 + 0.1210185*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.127005 + 0.1136015*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9711095 + 0.167852*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.590733 + 0.189764*x + -0.01111425*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.400373 + 0.137916*x + -0.009604315*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4380145 + 0.173627*x + -0.0108367*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4158195 + 0.127494*x + -0.008657875*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.382088 + 0.149924*x + -0.0100193*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.455317 + 0.2107825*x + -0.0133669*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5207225 + 0.193823*x + -0.01140495*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.3657465 + 0.198613*x + -0.0127062*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.521741 + 0.196093*x + -0.01131145*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.516351 + 0.160993*x + -0.01019485*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.416823 + 0.180059*x + -0.01064965*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(-0.646893 + 0.1334485*x + -0.0083625*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5936595 + 0.111935*x + -0.0080649*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(-0.629567 + 0.132426*x + -0.009647625*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.429283 + 0.1893055*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.222056 + 0.1935255*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2231245 + 0.1576465*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.1674005*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.255206 + 0.1621515*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.2759975 + 0.123709*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.374222 + 0.160055*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2011955*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(1.956205 + 0.142431*x + -0.004659135*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8367625 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4084025 + 0.1103895*x + -0.0101973*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.692884 + 0.138487*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9020235 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.521387 + 0.109307*x + -0.005736885*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
                stat_function(fun=function(x){(0 + 0.126315*x + -0.005174015*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.250244 + 0.194454*x + -0.006083915*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.3517415 + 0.20106*x + -0.00794801*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.378024 + 0.245699*x + -0.01145405*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1260425*x + -0.003445515*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4251285 + 0.07685265*x + -0.004730925*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.814752 + 0.13187*x + -0.008340265*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.773781 + 0.195075*x + -0.0122432*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.280487 + 0.125592*x + -0.006563*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.06244435*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.0630326*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(-0.705988 + 0.1707435*x + -0.007262235*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7657325 + 0.1303005*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6265235 + 0.1195185*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8095375 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9669255 + 0.1428735*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.07371 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-1.05278 + 0.1521105*x + -0.0045679*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.863773 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(-0.75 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5707205 + 0.2137465*x + -0.0106657*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.377741 + 0.261015*x + -0.0129277*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.815515 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.00993 + 0.0707026*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.897059 + 0.165041*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.002265 + 0.193276*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.639673 + 0*x + 0.009124395*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4271265 + 0.1881695*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.184442*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1861185*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.269866 + 0.270301*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.277012 + 0.3308325*x + -0.009847085*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.313874*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.215788*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.226289*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4703485 + 0.1919125*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8601045 + 0.125251*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.799565 + 0.149961*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7766935 + 0.1483915*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3990215 + 0.170607*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.372851 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.8730395 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3646365 + 0.116071*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.5780865 + 0.09902695*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.503424 + 0.08276695*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.415608 + 0.2575895*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.3535265 + 0.2517675*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.5701995 + 0.231955*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.7925675 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.807409 + 0.1946905*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.830252 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.647873 + 0.126091*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.947481 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6842465 + 0.102795*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.153812*x + -0.007501265*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1204285*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.139392*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.3389785 + 0.1869825*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.1735335*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.3587425 + 0.226134*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.939811 + 0.16052*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6155325 + 0.2556765*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.630523 + 0.1661375*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.617254 + 0.1522375*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.157309*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5125785 + 0.110644*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2205565*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.448908 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3208655 + 0.1520505*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2698865 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.69893 + 0.262698*x + -0.01650385*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.395351 + 0.2014525*x + -0.01070095*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.115518*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1468065*x + -0.01146465*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.2522065 + 0.09592255*x + -0.009500965*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1071145*x + -0.009444935*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.459635 + 0*x + -0.00956769*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.2538695 + 0.134041*x + -0.0112597*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.619206 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.727828 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3191835*x + -0.0109689*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.2791165*x + -0.0097306*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
                stat_function(fun=function(x){(0 + 0.2810145*x + -0.009506155*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
                              stat_function(fun=function(x){(0 + 0.3150845*x + -0.01127815*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(-0.838739 + 0.09976145*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.15033 + 0.1968865*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7507745 + 0.386003*x + -0.0106255*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.11676 + 0.178316*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-1.298205 + 0.08298065*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.12568 + 0.1978265*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.11508 + 0.132109*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.9910815 + 0.1844125*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-1.001215 + 0.1196895*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.9765795 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.013995 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.912906 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8685015 + 0.0993641*x + -0.008306955*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4604675 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-1.456075 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.67413 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.526694 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.541302 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.747209 + 0.107499*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.7339715 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4589845 + 0.1102405*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.579963 + 0.1262895*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.403334 + 0.1146835*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.533806 + 0.1857755*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.2524435 + 0*x + -0.0108022*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0107498*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.684806 + 0.1381225*x + -0.0146405*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6100035 + 0.09296875*x + -0.01490795*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.390845 + 0*x + -0.00866743*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-1.329735 + 0*x + -0.01008015*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.91243 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.004939005*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-1.649315 + 0.2199335*x + -0.01214475*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-0.53718 + 0*x + -0.008671945*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-1.076905 + 0.234262*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.083185 + 0.246818*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-1.0608 + 0.1684815*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.062665 + 0.2154275*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6313155 + 0.08813155*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.847431 + 0.147627*x + -0.00868069*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.753947 + 0.220392*x + -0.01124715*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8711005 + 0.144771*x + -0.008675095*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2215845*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.549339 + 0.1398475*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2732865*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.921255 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5786435 + 0.124246*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1095285*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4673615 + 0.182799*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3241165 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0108745*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.740998 + 0*x + -0.0129666*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.511347 + 0*x + -0.0117811*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.478088 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1687575*x + -0.01352165*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.01123455*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.199298*x + -0.0143934*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8887615 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3418295 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.622662 + 0.116858*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3404345 + 0.1198685*x + -0.00888389*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.3905855 + 0.14425*x + -0.009714605*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.424854 + 0.187779*x + -0.0115834*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
                stat_function(fun=function(x){(-0.23741 + 0.1921335*x + -0.01249275*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
                              stat_function(fun=function(x){(0 + 0.1654925*x + -0.01130075*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
                                            stat_function(fun=function(x){(0 + 0.172764*x + -0.0119523*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0.626957 + 0.184845*x + -0.01392835*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
                stat_function(fun=function(x){(0 + 0.139735*x + -0.0108491*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.2192755 + 0*x + 0*x^2)*(0.1701297)+(0.3140121)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.5441655 + 0.18185) + (0.132175+0.03107995)*x + 0*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.5441655 + 0.731963) + (0.132175 + 0.212783)*x + 0*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#EC1804')



# meanInsetPlot <- function() {
#   print(meanPlot)
#   
#   theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
#   
#   print(meanInset, vp=viewport(width = 0.13, height = 0.25, x = 0.115, y = 0.99, just = c("left","top")))
#   
#   theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
# }


print(meanPlot) #export at 1200x1000




###dispersion
chainsDispersionIntercept <- chainsCommunity[,297:586]%>%
  gather(key=parameter, value=value, B.1.2.1:B.290.2.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsDispersionSlope <- chainsCommunity[,1457:1746]%>%
  gather(key=parameter, value=value, B.1.2.2:B.290.2.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsDispersionQuad <- chainsCommunity[,2617:2906]%>%
  gather(key=parameter, value=value, B.1.2.3:B.290.2.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

dispersion <- read.csv('dispersion_experiment_coefs.csv')%>%
  left_join(expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  #get estimates at 10 years
  mutate(yr10=Intercepts + 10*Slopes + (10^2)*Quads,
         yr20=Intercepts + 20*Slopes + (20^2)*Quads,
         final_year_estimate=Intercepts + alt_length*Slopes + (alt_length^2)*Quads)%>%
  mutate(curve1='stat_function(fun=function(x){',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2}, size=0.5, xlim=c(0,',
         curve5='), colour=#bdbdbd) +',
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, alt_length, curve5, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

# #temporary solution until we get the model output for each plot mani
# dispersion2 <- dispersion%>%
#   summarise(mean_intercept=mean(Intercepts), mean_linear=mean(Slopes), mean_quadratic=mean(Quads), alt_length=max(alt_length))%>%
#   mutate(curve1='stat_function(fun=function(x){',
#          curve2=' + ',
#          curve3='*x + ',
#          curve4='*x^2}, size=3, xlim=c(0,',
#          curve5='), colour=black)',
#          curve=paste(curve1, mean_intercept, curve2, mean_linear, curve3, mean_quadratic, curve4, alt_length, curve5, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

# #change in dispersion at 10 years vs. MAT
# ggplot(data=dispersion, aes(x=MAT, y=yr10)) +
#   geom_point() +
#   stat_smooth(method='gam', se=T, formula = y ~ s(x), color='black') +
#   xlab(expression(atop('Mean Annual Temperature ' ( degree~C)))) +
#   ylab('Change in Dispersion at 10 Years') +
#   scale_y_continuous(limits=c(-0.4,0.6), breaks=seq(-0.4,0.6,0.2)) +
#   scale_x_continuous(breaks=seq(-50,20,10)) +
#   theme(legend.position='none')

#dispersion density plot with all data points
dispersionInset <- ggplot(data=dispersion, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('Dispersion Change') +
  ylab('Density')


#main figure
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.5,0.5))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){0.010688636 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.015890051*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.001485225*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.02191456 + -0.014966281*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.019016194*x + 0.001742704*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.064169496 + -0.017841245*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.048682274 + -0.016542292*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.015905597*x + 0*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.001327016*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.01386412*x + 0.001440809*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.014368443 + -0*x + 0.001283393*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.013431967*x + 0.001598169*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.019460086*x + 0.001825246*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.037927391 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.022398925 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027245251 + -0.028956762*x + 0.004383552*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.048314479 + -0.025357947*x + 0.004179582*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.060453159 + -0.02624632*x + 0.004076132*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.03776382 + -0.023438616*x + 0.004032436*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.092039053 + -0.023328118*x + 0.004064221*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026244749 + -0.034279431*x + 0.004698383*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016939834 + -0.030599171*x + 0.00448058*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.018613652 + -0.029626316*x + 0.004444485*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025015276 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029885487 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.062337365 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.034006512 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.03492158 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046096515 + 0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000733447*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.00038454*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000382396*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000270136*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000231877*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.011773665 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.013380707*x + 0.002480388*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.013567029*x + 0.002453951*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.002277486*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.014920641*x + 0.00265144*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.044805192 + -0.016957948*x + 0.002673399*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.001999204*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.040278664 + -0*x + 0.002435196*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002587254*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.069755551 + -0*x + 0.002555515*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.00248293*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.041189925 + -0.015925992*x + 0.002650321*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.045045358 + -0.014105102*x + 0.002588786*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002345792*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.068277392 + -0*x + 0.002489942*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.046407853 + -0.016131033*x + 0.002659608*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.002196522*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002537997*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002596487*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.002330419*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002316246*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,13), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,13), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026248964 + -0*x + -0*x^2}, size=0.5, xlim=c(0,13), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.080184064 + 0*x + 8.87e-06*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.210704821 + 0*x + 0.000188387*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.007903859*x + 0.000451697*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.096746283 + -0.015661624*x + 0.000118136*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.007695788 + -0.011028361*x + 0.000580176*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.066514953 + 0.017938236*x + -0.000965136*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.036551434 + 0*x + -0*x^2}, size=0.5, xlim=c(0,20), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,20), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.034844456 + -0*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.139358515 + -0.026485897*x + 0.000712483*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000899866*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.043426773 + 0*x + -0.000621176*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.013887114 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027205729 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017149679 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.001921788*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002532903*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002379028*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.002399083*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.001823152*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.001535579*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.011750233 + -0.022274045*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024704815 + -0.024483534*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.012731337 + -0.023498986*x + 0*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.054257708 + -0*x + 0.001238737*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.108294996 + 0*x + -0^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.078678711 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.054222689 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.030289238 + 0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.010842507 + -0*x + 0.00011896*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024459709 + -0*x + 0.000121224*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.012896992 + -0*x + 0.000286503*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027250281 + 0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000249605*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.028448572 + -0*x + 0.000213771*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000206285*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.000180665*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.000139687*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.020072503 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.013722456 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.182109736 + 0.045581451*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.177389815 + 0.038850692*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.138861776 + 0.037070321*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.104224703 + 0.036439517*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.13290001 + 0.046077374*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.182133304 + 0.044615214*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.076156359 + 0.04159499*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.109072887 + 0.04023526*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.03682103 + 0.008410305*x + -0*x^2}, size=0.5, xlim=c(0,23), colour='#bdbdbd') +
  stat_function(fun=function(x){0.013009709 + -0*x + 0.001002301*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.016588748*x + 0.001951543*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.032498213 + -0*x + 0.001744399*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046579883 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.047868003 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.089384041 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.056999468 + -0.01674407*x + 0.000791166*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.066203766 + -0.015095089*x + 0.000525694*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.01504824 + -0.022216259*x + 0.001066483*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025387558 + -0.027398744*x + 0.001303421*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.050425753 + -0.017869209*x + 0.000788277*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.057926773 + -0.01793887*x + 0.000846591*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.069577326 + -0.010052017*x + 0.000399116*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.035896837 + -0.017541434*x + 0.000772584*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.042974859 + -0.00850663*x + 0.000494566*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.062234981 + -0*x + 2.19e-05*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.057485465 + -0*x + 6.87e-05*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.059752755 + -0*x + 0.000543005*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.11228537 + -0*x + 0.000672972*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.011181703 + -0*x + 0.000463698*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000349383*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027026976 + -0*x + 0.000501144*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.000361442*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.05246031 + -0.014709754*x + 0.000989837*x^2}, size=0.5, xlim=c(0,18), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,18), colour='#bdbdbd') +
  stat_function(fun=function(x){0.015034417 + 0*x + 0.000176039*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.00042538*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025436733 + 0.014442848*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){0.041540165 + -0*x + 0.000611532*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017629602 + 0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000795488*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000739344*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016881232 + 0*x + -0.001444389*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.100564279 + 0.029826643*x + -0.002676858*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.105979452 + 0.012976427*x + -0.00194226*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0.001309907*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.009765032*x + -0.002048706*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.028660728 + -0*x + -0.001676102*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.013515089*x + 0.001269964*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.047812815 + -0.022503923*x + 0.001805503*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.013199688*x + 0.001548151*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029164492 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.079037925 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.061455428 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024120603 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026408183 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017405209 + -0.023347426*x + 0.003466966*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.074649133 + -0*x + 0.003530536*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.00366252*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.005222487*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.01776057*x + 0.005440173*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016268422 + -0.025023555*x + 0.005799574*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.006245197*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.03734395 + -0.019054265*x + 0.003203527*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.035384281 + -0.025677156*x + 0.003472826*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046574671 + -0.016795194*x + 0.003191933*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.018002776*x + 0.003359084*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.015894644 + -0.020252193*x + 0.003361754*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0.016060148*x + 0.003299181*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046419622 + -0*x + 0.003028268*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.003038629*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009105691 + -0.01787492*x + 0.003426049*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.038857234 + -0*x + 0.002989843*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.088326659 + -0*x + 0.002945228*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.078624037 + -0*x + 0.003223161*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.079755068 + -0*x + 0.003121896*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.084559788 + -0.016722496*x + 0.003319893*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.116166667 + -0*x + 0.002916833*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026439546 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.014587352 + -0.017683839*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.050088656 + -0*x + 0.00300267*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + 0*x + 0.002377777*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.00208873*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.101463031 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.149057904 + 0.020007223*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.146392921 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.015438137*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.023209645*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.032252065 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.054032287 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.04480628 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.072273688 + -0*x + 0.001729982*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.030507362 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.048003655 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.073751832 + 0.020799331*x + -0.001305697*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.041356607 + 0.01246301*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.101425413 + 0.021188019*x + -0.001387468*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016473553 + 0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.01722922 + 0.01510089*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.020483536 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.062099073 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.061477606 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024493338 + 0*x + -0.001446184*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017074624 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.01017463 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009177981 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029481117 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.041444352 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029051003 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0.027479827*x + 0.002832464*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.028850903 + -0.024488338*x + 0.002643486*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000474647*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.084593553 + -0*x + 0*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.038785035 + -0.016243252*x + 0.000715714*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.030663272 + -0*x + 0.000254553*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046179229 + -0.007734769*x + 0.000133632*x^2}, size=0.5, xlim=c(0,16), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0.01217027*x + 0.000738937*x^2}, size=0.5, xlim=c(0,16), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.015141046*x + -0*x^2}, size=0.5, xlim=c(0,16), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + 0.000527392*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0 + -0*x + 0.00119602*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.042078283 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025795056 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.047775364 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.068824756 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.014821747*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.052934749 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.09633173 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.022549336*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.072243956 + 0.015875412*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + 0.020341071*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.00836217 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.032822921 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.037871915 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.008411798 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.02892351 + 0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.058116585 + 0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.020997814 + 0.015313182*x + -0.002257762*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009484363 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.044456841 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.039234879 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.066933796 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.056481307 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.056221637 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.104099159 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.053806444 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.070187613 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.031836741 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.060027139 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.126135123 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.079409497 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.227560838 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.131284251 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.064842918 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  #estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0 - 0}, size=3, xlim=c(0,23), colour='black')


dispersionInsetPlot <- function() {
  print(dispersionPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(dispersionInset, vp=viewport(width = 0.13, height = 0.25, x = 0.83, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

dispersionInsetPlot() #export at 1200x1000







###richness
chainsRichnessIntercept <- chainsCommunity[,877:1166]%>%
  gather(key=parameter, value=value, B.1.4.1:B.290.4.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsRichnessSlope <- chainsCommunity[,2037:2326]%>%
  gather(key=parameter, value=value, B.1.4.2:B.290.4.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsRichnessQuad <- chainsCommunity[,3197:3486]%>%
  gather(key=parameter, value=value, B.1.4.3:B.290.4.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)


richness <- read.csv('richness_experiment_coefs.csv')%>%
  left_join(expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  #get estimates at 10 years
  mutate(yr10=Intercepts + 10*Slopes + (10^2)*Quads,
         yr20=Intercepts + 20*Slopes + (20^2)*Quads,
         final_year_estimate=Intercepts + alt_length*Slopes + (alt_length^2)*Quads)%>%
  mutate(curve1='stat_function(fun=function(x){',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

# #temporary solution until we get the model output for each plot mani
# richness2 <- richness%>%
#   group_by(plot_mani)%>%
#   summarise(mean_intercept=mean(Intercepts), mean_linear=mean(Slopes), mean_quadratic=mean(Quads), alt_length=max(alt_length))%>%
#   mutate(curve1='stat_function(fun=function(x){',
#          curve2=' + ',
#          curve3='*x + ',
#          curve4='*x^2}, size=3, xlim=c(0,',
#          curve5='), colour=',
#          curve6=') +',
#          color=ifelse(plot_mani==1, '#1400E5', ifelse(plot_mani==2, '#4A06AC', ifelse(plot_mani==3, '#800C74', ifelse(plot_mani==4, '#B6123C', '#EC1804')))),
#          curve=paste(curve1, mean_intercept, curve2, mean_linear, curve3, mean_quadratic, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

# #change in richness at 10 years vs. MAT
# ggplot(data=richness, aes(x=MAT, y=yr10)) +
#   geom_point() +
#   # stat_smooth(method='gam', se=T, formula = y ~ s(x), color='black') +
#   xlab(expression(atop('Mean Annual Temperature ' ( degree~C)))) +
#   ylab('Proportional Change in Richness at 10 Years') +
#   scale_y_continuous(limits=c(-1.5,1.5), breaks=seq(-2,2,0.5)) +
#   scale_x_continuous(breaks=seq(-50,20,10)) +
#   theme(legend.position='none')



#inset
richnessInset <- ggplot(data=richness, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('Proportion Richness Change') +
  ylab('Density')


#main figure
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-1,1))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.018172345 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.033655242*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.003885528*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.010734327 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.07652313 + -0*x + -0.000376624*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.081341003 + 0*x + -0.002606891*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.062401703 + -0*x + -0.001075759*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002493991*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.159103355 + -0*x + 0.003843*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.20294003 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.169927674 + -0.029430969*x + -0.004460169*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.190826393 + -0*x + -0.003500128*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.20369155 + -0*x + 0.011081314*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.046223084*x + 0.00793682*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.00848739*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.005802431*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.008659413 + -0.039544248*x + 0.007704011*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.149098024 + -0*x + 0.01058593*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.161159402 + -0*x + 0.010050752*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.192147157 + -0*x + 0.010303687*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.041279427 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.057459101 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011299454 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.065861648 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.02295008 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.033512305 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.042908483 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.045069962 + -0.031694564*x + 0.001287275*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.003253768*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.048015667 + -0.027549195*x + 0.001426462*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.136630589 + -0*x + 0*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.020796263 + 0.027382012*x + -0.005741985*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0.04016861*x + -0.005971992*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0*x + -0.004118268*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0.028087563*x + -0.005538724*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0.031221262*x + -0.00556602*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.046603421 + 0.053060982*x + -0.006878391*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.07964425 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.009718768 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.033597494*x + 0.003569036*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.032316639*x + 0.002962308*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002761689*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.002887236*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.019965839 + 0*x + -0.000755811*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.00318964*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002872976*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.09698491 + -0.028272585*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.122271653 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.051309 + -0.036925819*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.004237366*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.062389125 + -0.027518205*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.0205894 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.042744042*x + 0.00375886*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.004263495*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.028799861*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.019964087*x + 0.001978936*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.002516999*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){-0.328470365 + -0*x + 0.004448551*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.12349762 + 0*x + -0.000850106*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002403951*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.031874443 + -0*x + 0*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){-0.282024181 + -0*x + 0.003265135*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){-0.286311118 + -0*x + 0.002679653*x^2}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){-0 + -0*x + 0.003486886*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){-0.331620759 + -0*x + 0.00705288*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.013641266 + -0.021222563*x + 0.00193886*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){-0.209323951 + -0*x + 0.006384082*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){-0 + -0.022822678*x + 0.002679459*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){-0.305382265 + -0*x + 0.005974897*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.09024787 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.06409514 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.028767503 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.034525119 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.044584879*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.044974133*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.005701355*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.196543078 + -0*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.088434821 + -0*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.038273002 + -0.034161592*x + 0*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0*x + 0.003903584*x^2}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){-0.21376675 + -0*x + 0.003280229*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.188736841 + -0*x + 0.003020113*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.062438067 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.0042668*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.122903448 + -0.033134941*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.03803995*x + 0.004555676*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.026624549 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.126663337 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.009614823 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.042036322 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011881277 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.089913963 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0.02623026*x + -0.002507142*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0.001913016*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.01397403 + 0*x + -0.001230633*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.070134863 + 0.020158885*x + -0.00261828*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.037242463 + 0*x + -0.001823333*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0.018707886*x + -0.002298656*x^2}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + 0.024146655*x + -0.002120516*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + -0.001407808*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.033891953 + 0*x + -0.002176754*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.027543709 + 0.022690327*x + -0.002711683*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.009819284 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0.02868894*x + -0.002484557*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0.0015047*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.055875871 + 0.041780628*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.094887811 + 0.033769593*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.043857604 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.023367579 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.074068224 + 0.032359749*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.080227902 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.045932633 + -0*x + -0*x^2}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){-0.130744892 + -0*x + 0.00168356*x^2}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.031452253*x + 0.004394314*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.007003369 + -0*x + 0.005546626*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.046889298*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.003076305*x^2}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002846664*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0*x + 0.002755388*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0*x + 0.005854724*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0.145142822 + -0*x + 0.004016289*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002592162*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.025218704*x + 0.00122207*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.026503676*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.031110912*x + -0.002517045*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.176030405 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.117342402 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011234731 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.02085897*x + 0*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){0.090323114 + 0.033779885*x + -0.003453232*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.094726918 + 0.01882863*x + -0.002339971*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.013808106 + -0*x + -0.001996972*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.032952435*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.043315462 + 0.025278853*x + -0.004092482*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.028524731*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.023703545 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.031836132*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.144074095 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.160032813 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.124682191 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.160997943 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.089702183 + -0.054107874*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.145877767 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.069240288 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.152396504 + -0*x + 0.007380956*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0.040973075*x + 0.007602615*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.041330082*x + 0.008157942*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.234994002 + -0*x + 0.008411243*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0.181823928 + -0*x + 0.008017049*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.173282301 + -0.039074491*x + 0.007135362*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.037526167*x + 0.005962066*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.009246229*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.056066839 + -0*x + 0.007183138*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.026831408*x + 0.007177524*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.007197329*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0.006479669*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.044464238*x + 0.007257364*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.008681354*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.035463948 + -0*x + 0.008942751*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.039765294 + -0*x + 0.0080208*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.009297573*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.046604104*x + 0.008202677*x^2}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0.028987741*x + 0.006762392*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.01897864 + -0.045453839*x + 0.007324294*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.027172866 + -0.036390653*x + 0.006757898*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.164054568 + -0.032235677*x + 0.005631738*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.293921699 + 0.098790523*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.119439043 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.292126374 + 0.103761056*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.155884807 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.0439309*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.018066684 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.023848428 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.082439353 + -0.0378053*x + -0.003712491*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.160746369 + 0*x + -0.006900301*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.042242724*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.038721487*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.038785752*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.076284135 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.041465245*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.027044217*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0*x + 0.00247231*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0*x + 0.002385414*x^2}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.028350461 + 0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.156562234 + 0*x + -0.003569162*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.044519015 + 0.024257348*x + -0.002740088*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.017182286 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.029049151 + -0.03147215*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.11045396 + 0*x + -0.002493054*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.028863212*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.063627327 + 0*x + -0.003590923*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.031505773 + 0.043652454*x + -0.004862013*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){0.088315667 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.097580112 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.058520258 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.075802942 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.093180424 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.05075954 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.035191333 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.139736951 + -0*x + 0.006835795*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.026462969 + -0.034745815*x + 0.004697403*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0.003310635*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.168289018 + 0.029145776*x + -0.002974509*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.037702383*x + 0.001956475*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.116965155 + 0*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){0.050216597 + 0*x + -0*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.057770567 + 0*x + -0*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.041115786*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.078386837 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0.031212343*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.032259923 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.059198022 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.101819479 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.096039907 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.116357375 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.027886418*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0.028408206*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.024395488 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.079791421 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.068256748 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.115979442 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.073630943 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.083969266 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.03844796 + -0.039719321*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.197549374 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.009026597 + -0.034835696*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.046179058*x + -0.005080928*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.110014848 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.041992639 + 0.034363994*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.028205006 + 0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.019189047 + -0*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.042602249 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.007232073 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){-0.166418787 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){-0.300856218 + -0*x + 0.003821639*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
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
                                  


richnessInsetPlot <- function() {
  print(richnessPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(richnessInset, vp=viewport(width = 0.13, height = 0.25, x = 0.30, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

richnessInsetPlot() #export at 1200x1000






###evenness
chainsEvennessIntercept <- chainsCommunity[,587:876]%>%
  gather(key=parameter, value=value, B.1.3.1:B.290.3.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsEvennessSlope <- chainsCommunity[,1747:2036]%>%
  gather(key=parameter, value=value, B.1.3.2:B.290.3.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsEvennessQuad <- chainsCommunity[,2907:3196]%>%
  gather(key=parameter, value=value, B.1.3.3:B.290.3.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)


evenness <- read.csv('evenness_experiment_coefs.csv')%>%
  left_join(expInfo, by=c('site_code', 'project_name', 'community_type', 'treatment'))%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  #get estimates at 10 years
  mutate(yr10=Intercepts + 10*Slopes + (10^2)*Quads,
         yr20=Intercepts + 20*Slopes + (20^2)*Quads,
         final_year_estimate=Intercepts + alt_length*Slopes + (alt_length^2)*Quads)%>%
  mutate(curve1='stat_function(fun=function(x){',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2}, size=0.5, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E544', ifelse(plot_mani==2, '#4A06AC44', ifelse(plot_mani==3, '#800C7444', ifelse(plot_mani==4, '#B6123C44', '#EC180444')))),
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

# #temporary solution until we get the model output for each plot mani
# evenness2 <- evenness%>%
#   group_by(plot_mani)%>%
#   summarise(mean_intercept=mean(Intercepts), mean_linear=mean(Slopes), mean_quadratic=mean(Quads), alt_length=max(alt_length))%>%
#   mutate(curve1='stat_function(fun=function(x){',
#          curve2=' + ',
#          curve3='*x + ',
#          curve4='*x^2}, size=3, xlim=c(0,',
#          curve5='), colour=',
#          curve6=') +',
#          color=ifelse(plot_mani==1, '#1400E5', ifelse(plot_mani==2, '#4A06AC', ifelse(plot_mani==3, '#800C74', ifelse(plot_mani==4, '#B6123C', '#EC1804')))),
#          curve=paste(curve1, mean_intercept, curve2, mean_linear, curve3, mean_quadratic, curve4, alt_length, curve5, color, curve6, sep='')) #need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below

# #change in evenness at 10 years vs. gamma diversity
# ggplot(data=evenness, aes(x=rrich, y=yr10)) +
#   geom_point() +
#   # stat_smooth(method='gam', se=T, formula = y ~ s(x), color='black') +
#   xlab('Gamma Diversity') +
#   ylab('Change in Evenness at 10 Years') +
#   scale_y_continuous(limits=c(-0.3,0.6), breaks=seq(-0.4,0.6,0.1)) +
#   scale_x_continuous(breaks=seq(0,200,25)) +
#   theme(legend.position='none')




#inset
evennessInset <- ggplot(data=evenness, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('Evenness Change') +
  ylab('Density')


#main figure
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.4,0.6))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){0 + 0.02967549*x + -0.001713149*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0.025216304*x + -0.001405815*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.086278093 + 0.017750593*x + -0.001452646*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.082879172 + 0*x + -0.001448261*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.05019695 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.012809289*x + 0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0.043284622 + -0.031154267*x + 0.001396124*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.000897488*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.015086716*x + -0.002430533*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.013395711*x + -0.002329999*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + -0.002120862*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008439051 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.016833005 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.0112123 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.046155757 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.047001726 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.009044503 + -0.021786924*x + 0.001164979*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){-0.064355067 + -0.031615234*x + 0.001953928*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.028336341 + -0.011992486*x + 0.00066197*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.012840988*x + 0.002229226*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0.015043846*x + 0.002526124*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0.001765482*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0.00209756*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.022541209 + -0.025946022*x + 0.003086287*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0.013275375*x + 0.002288334*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.038968628 + 0.011550912*x + -0.001770808*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.049013614 + 0*x + -0.001458624*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.072974282 + 0*x + -0.001155096*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.049905348 + 0.011499509*x + -0.001794348*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.067936379 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001330982*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001163965*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.00132096*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001265012*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0.001271145*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001363592*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001376054*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001175952*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001411509*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.010766792*x + -0.001445898*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.013417637*x + -0.001848129*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.017877861 + 0*x + -0.001596785*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.017021237*x + -0.002126633*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.018973598*x + -0.002400245*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.014169956 + 0*x + -0.001535764*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.022042385*x + -0.00170563*x^2}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.021564663*x + -0.001696978*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){-0.018244659 + 0.012914834*x + -0.000809716*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.011239187*x + -0.00043966*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.070675543 + 0.05034904*x + -0.003066594*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){-0.014035141 + -0*x + 0*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.035473984*x + -0.001257482*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.032369929*x + -0.001469313*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.047965237 + 0*x + -0.000463148*x^2}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.017518661*x + -0.000870656*x^2}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.164597191 + 0.093763562*x + -0.006221038*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){-0.021927586 + 0*x + -0.00068973*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.065170315*x + -0.003456405*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0 + 0.008147414*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.048082421*x + -0.00239693*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.030067062 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.025454444 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){-0.007720244 + 0.016154945*x + -0.00180316*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.020472232*x + -0.002132239*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.013916616 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.009266558 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.0200077 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0.001321803*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0*x + -0.001321787*x^2}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){0.046735871 + 0.02387679*x + -0.002059546*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.013622885*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.029867192 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.035814619 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.03645022 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.062646686 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.060674105 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.05696377 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.056744347 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.047267167 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.010403562*x + 0.000516358*x^2}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.012132111*x + 0.00046048*x^2}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.011040611*x + -0.001018194*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001069481*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.048183884 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.046039255 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.000546895*x^2}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){0.049862421 + -0*x + -0*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.048795743 + -0*x + -0*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0.008299449*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0.011563426*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.048663255 + -0*x + -0*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.008048059*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.013655774*x + 0.000619645*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.014781456*x + 0.001129917*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.018631622*x + -0.001071606*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.016244411 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.018130836*x + -0.001278183*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.019056649*x + -0.001132343*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.014033121*x + -0.000770548*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.006828601*x + 0*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.012033701*x + -0.000747588*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.045413215 + -0*x + 0*x^2}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.056572351 + 0.018205406*x + -0.001539311*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.015153584 + 0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.028610212 + 0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.01262416 + 0*x + -7.6e-05*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.057744264 + -0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.034889892 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.04448384 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.030402005 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){-0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.015276039 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.020789182 + -0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0.010672936 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.020965833 + -0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.018556173*x + -0.001466841*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + 0.012599904*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.016012945 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0.001296991*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0*x + -0.001293919*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001125362*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.013392556*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.067314062 + -0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0.014361756*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.011229727 + -0*x + 0*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.12145939 + 0.020537545*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.118716031 + 0.014361601*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.115461233 + 0.013045564*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.021438277*x + 0.001422074*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.020463286*x + 0.001381086*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.015558649 + 0.030534594*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.054101107*x + -0.002664611*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.058818337 + 0.059377975*x + -0.003073234*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.060921803 + 0.089161183*x + -0.004407644*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.009822883 + 0.03309308*x + -0.001419524*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.037871127*x + -0.001585837*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.012493586 + 0.036455543*x + -0.001489279*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008640162 + 0.03457108*x + -0.001325541*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.037705698*x + -0.001405996*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.04268077*x + -0.001782785*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.044984623*x + -0.001941207*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.016010719*x + -0.001256386*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.023346373*x + -0.001558192*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.024587064*x + -0.001550011*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.025985582*x + -0.00170006*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.021055758*x + -0.001329353*x^2}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){0 + 0.016266231*x + -0.000981826*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + 0.001635781*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.048633988*x + -0.002858259*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0.00170141*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.027548495*x + 0.002674845*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.03411544*x + 0.003109371*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.021543926*x + 0.002182266*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.009431107 + -0.033002789*x + 0.003116068*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){-0.022996074 + -0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008113209 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.007312179 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.051524835 + -0*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.049843448 + -0*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.043106886 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){-0.030977878 + -0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){-0.012478175 + 0.013712658*x + -0.001519867*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.018334326*x + -0.001808564*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.011853896*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.013938435*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.036546302 + 0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.022600678 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.013712348*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){-0.012468505 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008068891 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.010823055*x + -0.000884698*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0. + 0.022532363*x + -0.001136707*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){-0.032407889 + -0*x + 0*x^2}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.009019021*x + -0*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.07658946 + 0.021887919*x + -0.00146121*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0.010367987 + 0.011320014*x + -0.000805368*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.009392951*x + -0.00081311*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.015245477 + 0.017808097*x + -0.001890435*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.023020482 + 0.012028427*x + -0.001438751*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0 + 0.021117658*x + -0.002229459*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.011432751 + 0.021566576*x + -0.002105576*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.027393205 + 0.021872871*x + -0.002174845*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.033500217 + 0.012077554*x + -0.001548165*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.022029305 + 0.015568523*x + -0.001991342*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.053827705 + 0*x + -0.001431278*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.031448826 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.023485325*x + 0.001976154*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + -0.032044553*x + 0.002628645*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.04645807 + -0.027919844*x + 0.001982293*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0.029427476*x + 0.002424643*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.020519863 + -0*x + 0.001389676*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.035720915 + 0*x + -0.001389324*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.023960926*x + -0.002109459*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){-0.015995363 + 0.021007355*x + -0.001828975*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.044195593 + 0.021443759*x + -0.0017091*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.036319542 + 0.052472105*x + -0.003171142*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.012459455 + 0.025667652*x + -0.002161746*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){-0.042320336 + 0.048739534*x + -0.00292485*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.046939619 + 0.021907105*x + -0.002626824*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.046464844 + 0.016370484*x + -0.002491818*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001828142*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.001676888*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.048482777 + 0.017583871*x + -0.002622398*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.047044802 + 0.017439179*x + -0.002559794*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.062033316 + 0.025574466*x + -0.003014897*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + 0*x + -0.001580578*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.052528883 + 0.019032602*x + -0.002512958*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0 + 0*x + -0.001980907*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0 + 0*x + -0.001339912*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0 + 0*x + -0.001943581*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
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



evennessInsetPlot <- function() {
  print(evennessPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(evennessInset, vp=viewport(width = 0.13, height = 0.25, x = 0.80, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

evennessInsetPlot() #export at 1200x1000



# ###all four variables in one figure
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanInsetPlot(), vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionInsetPlot(), vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessInsetPlot(), vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessInsetPlot(), vp=viewport(layout.pos.row = 2, layout.pos.col = 2))







###magnitude change in the four variables at 10 yrs, 20 yrs, final year
meanEstimates <- mean%>%
  summarise(yr10_median=median(yr10), yr20_median=median(yr20), final_median=median(final_year_estimate), yr10_sd=sd(yr10), yr20_sd=sd(yr20), final_sd=sd(final_year_estimate))%>%
  mutate(yr10_CI=yr10_sd*2, yr20_CI=yr20_sd*2, final_CI=final_sd*2)%>%
  gather(key=timepoint, value=estimate)%>%
  separate(timepoint, c('timepoint', 'stat'), sep = '_')%>%
  spread(key=stat, value=estimate)%>%
  mutate(variable='mean change')
dispersionEstimates <- dispersion%>%
  summarise(yr10_median=median(yr10), yr20_median=median(yr20), final_median=median(final_year_estimate), yr10_sd=sd(yr10), yr20_sd=sd(yr20), final_sd=sd(final_year_estimate))%>%
  mutate(yr10_CI=yr10_sd*2, yr20_CI=yr20_sd*2, final_CI=final_sd*2)%>%
  gather(key=timepoint, value=estimate)%>%
  separate(timepoint, c('timepoint', 'stat'), sep = '_')%>%
  spread(key=stat, value=estimate)%>%
  mutate(variable='dispersion')
richnessEstimates <- richness%>%
  summarise(yr10_median=median(yr10), yr20_median=median(yr20), final_median=median(final_year_estimate), yr10_sd=sd(yr10), yr20_sd=sd(yr20), final_sd=sd(final_year_estimate))%>%
  mutate(yr10_CI=yr10_sd*2, yr20_CI=yr20_sd*2, final_CI=final_sd*2)%>%
  gather(key=timepoint, value=estimate)%>%
  separate(timepoint, c('timepoint', 'stat'), sep = '_')%>%
  spread(key=stat, value=estimate)%>%
  mutate(variable='richness')
evennessEstimates <- evenness%>%
  summarise(yr10_median=median(yr10), yr20_median=median(yr20), final_median=median(final_year_estimate), yr10_sd=sd(yr10), yr20_sd=sd(yr20), final_sd=sd(final_year_estimate))%>%
  mutate(yr10_CI=yr10_sd*2, yr20_CI=yr20_sd*2, final_CI=final_sd*2)%>%
  gather(key=timepoint, value=estimate)%>%
  separate(timepoint, c('timepoint', 'stat'), sep = '_')%>%
  spread(key=stat, value=estimate)%>%
  mutate(variable='evenness')
estimates <- rbind(meanEstimates, dispersionEstimates, richnessEstimates, evennessEstimates)
  
estimates10yr <- ggplot(data=subset(estimates, timepoint=='yr10'), aes(x=variable, y=median)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  ylab('Response at 10 years') +
  scale_x_discrete(limits=c('mean change', 'dispersion', 'richness', 'evenness'),
                   labels=c('Mean Change', 'Dispersion Change', 'Proportion Richness Change', 'Evenness Change')) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(angle=90, hjust=1))

estimates20yr <- ggplot(data=subset(estimates, timepoint=='yr20'), aes(x=variable, y=median)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  ylab('Response at 20 years') +
  scale_x_discrete(limits=c('mean change', 'dispersion', 'richness', 'evenness'),
                   labels=c('Mean Change', 'Dispersion Change', 'Proportion Richness Change', 'Evenness Change')) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(angle=90, hjust=1))

estimatesFinalyr <- ggplot(data=subset(estimates, timepoint=='final'), aes(x=variable, y=median)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  ylab('Response in Final Year') +
  scale_x_discrete(limits=c('mean change', 'dispersion', 'richness', 'evenness'),
                   labels=c('Mean Change', 'Dispersion Change', 'Proportion Richness Change', 'Evenness Change')) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(angle=90, hjust=1))
  
pushViewport(viewport(layout=grid.layout(1,3)))
print(estimates10yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(estimates20yr, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(estimatesFinalyr, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))


#final year estimates boxplot
meanEstimates <- mean%>%
  select(Intercepts, Quads, Slopes, final_year_estimate)%>%
  mutate(variable='mean change')
dispersionEstimates <- dispersion%>%
  select(Intercepts, Quads, Slopes, final_year_estimate)%>%
  mutate(variable='dispersion')
richnessEstimates <- richness%>%
  select(Intercepts, Quads, Slopes, final_year_estimate)%>%
  mutate(variable='richness')
evennessEstimates <- evenness%>%
  select(Intercepts, Quads, Slopes, final_year_estimate)%>%
  mutate(variable='evenness')
estimatesAll <- rbind(meanEstimates, dispersionEstimates, evennessEstimates, richnessEstimates)

ggplot(data=estimatesAll, aes(x=variable, y=final_year_estimate)) +
  geom_boxplot() +
  ylab('Response in Final Year') +
  scale_x_discrete(limits=c('mean change', 'dispersion', 'richness', 'evenness'),
                   labels=c('Mean Change', 'Dispersion Change', 'Proportion Richness Change', 'Evenness Change')) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(angle=45, hjust=1)) +
  geom_hline(aes(yintercept=0)) +
  annotate('text', x=1.1, y=0.46, label='*', size=10, hjust='left')
#export at 800x1000



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
