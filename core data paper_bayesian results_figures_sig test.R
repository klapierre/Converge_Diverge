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
#gather the intercepts, linear slopes, and quadratic slopes for all treatments
#set any that are not significant (CI overlaps 0) as 0
chainsDispersionIntercept <- chainsCommunity[,297:586]%>%
  gather(key=parameter, value=value, B.1.2.1:B.290.2.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-sd, upper=intercept+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsDispersionIntercept)[1] <- 'parameter1'

chainsDispersionSlope <- chainsCommunity[,1457:1746]%>%
  gather(key=parameter, value=value, B.1.2.2:B.290.2.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-sd, upper=slope+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsDispersionSlope)[1] <- 'parameter2'

chainsDispersionQuad <- chainsCommunity[,2617:2906]%>%
  gather(key=parameter, value=value, B.1.2.3:B.290.2.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-sd, upper=quad+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
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
  mutate(yr10=intercept + 10*slope + (10^2)*quad,
         yr20=intercept + 20*slope + (20^2)*quad,
         final_year_estimate=intercept + alt_length*slope + (alt_length^2)*quad)%>%
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


# #dispersion density plot with all data points
# dispersionInset <- ggplot(data=dispersion, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Dispersion Change') +
#   ylab('Density')


#main figure
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.35,0.35))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1752985*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.01638495*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1651075*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.209786*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.733904 + -0.196824*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5630495 + -0.182494*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.17547*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0146396*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1529485*x + 0.01589495*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.01415835*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.148181*x + 0.01763095*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.214683*x + 0.02013605*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.444402 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.31945*x + 0.0483592*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.558992 + -0.279748*x + 0.046109*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.6929055 + -0.2895485*x + 0.04496775*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.4425975 + -0.258574*x + 0.0444857*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(1.04136 + -0.257355*x + 0.04483635*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3781695*x + 0.0518324*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.337569*x + 0.0494296*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3268365*x + 0.0490314*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.355684 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.713692 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.401147 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.411242 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.5345235 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.008091365*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
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
  stat_function(fun=function(x){(0 + -0.1476155*x + 0.02736355*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.149671*x + 0.0270719*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02512515*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.164604*x + 0.0292506*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.468301 + -0.1870795*x + 0.02949285*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02205515*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4183645 + 0*x + 0.026865*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0285425*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7435525 + 0*x + 0.02819235*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0273916*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4284175 + -0.175695*x + 0.02923825*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4709505 + -0.155607*x + 0.0285594*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0258787*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.7272455 + 0*x + 0.02746895*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.4859815 + -0.177957*x + 0.0293407*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02423195*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0279991*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02864435*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0257091*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02555275*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,13), colour='grey') +
  stat_function(fun=function(x){(-0.8585995 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(-2.2985 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + -0.0871951*x + 0.00498311*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(1.09329 + -0.1727785*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1216645*x + 0.006400485*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0.759779 + 0.197894*x + -0.01064735*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0.4292225 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,20), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,20), colour='grey') +
  stat_function(fun=function(x){(-0.3584145 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(-1.51141 + -0.2921915*x + 0.00786009*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + -0.390211*x + 0.009927295*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.505071 + 0*x + -0.00685279*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0279429*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02624535*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0264666*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02011295*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2457265*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2701015*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.25924*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0.6245575 + -0.491978*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(1.220695 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.8419925 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.5721945 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.3081615 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.295827 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.3266125 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
  stat_function(fun=function(x){(0.339832 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,12), colour='grey') +
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
  stat_function(fun=function(x){(-0.38022 + 0.0927822*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,23), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3100105*x + 0.01105735*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,22), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1830065*x + 0.02152935*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.3845075 + 0*x + 0.01924415*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.539856 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.5540665 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(1.01207 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.6548045 + -0.18472*x + 0.008728115*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.756346 + -0.1665285*x + 0.00579944*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.245089*x + 0.0117654*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0 + -0.302262*x + 0.0143793*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.5822835 + -0.1971325*x + 0.00869624*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.6650345 + -0.197901*x + 0.00933956*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.793563 + -0.1108935*x + 0.00440303*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.422001 + -0.1935165*x + 0.008523125*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.5000855 + -0.09384485*x + 0.00545604*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.7125625 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(0.660166 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,21), colour='grey') +
  stat_function(fun=function(x){(-0.633202 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.21274 + 0*x + 0.0074242*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.552752 + -0.1622775*x + 0.01091985*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,18), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0 + 0.159333*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,11), colour='grey') +
  stat_function(fun=function(x){(0.484258 + 0*x + 0.006746405*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.008775795*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.00815642*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.01593445*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.13541 + 0.3290465*x + -0.029531*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(1.19515 + 0.1431555*x + -0.02142695*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.01445085*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1077275*x + -0.02260125*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0.3421725 + 0*x + -0.0184907*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + -0.149098*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.501481 + -0.2482625*x + 0.01991825*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1456185*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
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
  stat_function(fun=function(x){(0 + -0.195934*x + 0.0600158*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.276059*x + 0.0639807*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.4161715*x + 0.0688968*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.4379655 + -0.210206*x + 0.0353412*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.4163465 + -0.2832695*x + 0.0383121*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0.5397985 + -0.185284*x + 0.0352133*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.198606*x + 0.0370573*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2234215*x + 0.03708675*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + -0.177175*x + 0.03639645*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.538088 + 0*x + 0.03340775*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.03352205*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1971955*x + 0.03779605*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.45466 + 0*x + 0.03298385*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.000405 + 0*x + 0.03249165*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.893366 + 0*x + 0.0355578*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9058435 + 0*x + 0.03444065*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.958849 + -0.184482*x + 0.03662495*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.307535 + 0*x + 0.0321784*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + -0.1950875*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.526588 + 0*x + 0.03312535*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.02623155*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.0230428*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.145325 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.67039 + 0.220719*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(1.64099 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.170313*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0.256048*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.381792 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(-0.570094 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.468313 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7713325 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.3625445 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.555563 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.839616 + 0.2294575*x + -0.0144044*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0.482233 + 0.1374915*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(1.14491 + 0.2337455*x + -0.0153065*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1665925*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,10), colour='grey') +
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
  stat_function(fun=function(x){(0 + 0*x + -0.01595425*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,7), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.351223 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.483201 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.346478 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + -0.3031565*x + 0.03124765*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,5), colour='grey') +
  stat_function(fun=function(x){(0 + -0.2701545*x + 0.02916285*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,3), colour='grey') +
  stat_function(fun=function(x){(0.9592215 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.4538635 + -0.179195*x + 0.00789573*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.3642645 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,9), colour='grey') +
  stat_function(fun=function(x){(0.535436 + -0.0853297*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + -0.134262*x + 0.008151925*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1670355*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,16), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0.005818165*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,8), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.4901945 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0.5530445 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.733284 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(0 + 0.163513*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.557986 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.03674 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.2487635*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.7710045 + 0.175137*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0.224402*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.336113 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.44379 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.3450715 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0.6671285 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0.1689345*x + -0.02490755*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + -0.02288645*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.464458 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.4068495 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.712423 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(-0.5971115 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,4), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0.6462235 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.12243 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.5676025 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,6), colour='grey') +
  stat_function(fun=function(x){(-0.748319 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.3252335 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.636229 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.36553 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.8500545 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-2.484455 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-1.422335 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  stat_function(fun=function(x){(-0.6893565 + 0*x + 0*x^2)*(0.09064568)+(-0.00235573)}, size=0.5, xlim=c(0,2), colour='grey') +
  
  #estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*0 - 0}, size=3, xlim=c(0,23), colour='black')


# dispersionInsetPlot <- function() {
#   print(dispersionPlot)
#   
#   theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
#   
#   print(dispersionInset, vp=viewport(width = 0.13, height = 0.25, x = 0.83, y = 0.99, just = c("left","top")))
#   
#   theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
# }

print(dispersionPlot) #export at 1200x1000







###richness
chainsRichnessIntercept <- chainsCommunity[,877:1166]%>%
  gather(key=parameter, value=value, B.1.4.1:B.290.4.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-sd, upper=intercept+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsRichnessIntercept)[1] <- 'parameter1'

chainsRichnessSlope <- chainsCommunity[,2037:2326]%>%
  gather(key=parameter, value=value, B.1.4.2:B.290.4.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-sd, upper=slope+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsRichnessSlope)[1] <- 'parameter2'

chainsRichnessQuad <- chainsCommunity[,3197:3486]%>%
  gather(key=parameter, value=value, B.1.4.3:B.290.4.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-sd, upper=quad+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
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
  mutate(yr10=intercept + 10*slope + (10^2)*quad,
         yr20=intercept + 20*slope + (20^2)*quad,
         final_year_estimate=intercept + alt_length*slope + (alt_length^2)*quad)%>%
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



# #inset
# richnessInset <- ggplot(data=richness, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Proportion Richness Change') +
#   ylab('Density')


#main figure
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-1,0.65))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Proportion Richness Change')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4186895 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1471565*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.01698935*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.386167 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.340586 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6738265 + -0.143329*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0.6948925 + 0*x + -0.01139855*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0.2386625 + -0.2268595*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.612081 + -0.0873918*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.297744*x + 0.0109049*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.356443 + -0.3895715*x + 0.0168034*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5481175 + -0.2827385*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.403772 + -0.128686*x + -0.01950195*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.495151 + -0.180448*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5514035 + -0.3887175*x + 0.0484527*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.202109*x + 0.0347035*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2274025*x + 0.03711085*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.2831275 + 0*x + 0.02537095*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.3770945 + -0.172906*x + 0.03368555*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.312695 + -0.336831*x + 0.04628665*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.365433 + -0.319147*x + 0.0439466*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.500926 + -0.3267065*x + 0.04505255*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.5197245 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.388638 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6272095 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.3321525 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.43958 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.485763 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5268475 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.5362985 + -0.1385835*x + 0.00562857*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0.294763 + -0.239153*x + 0.014227*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5491785 + -0.120458*x + 0.00623716*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.9366445 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0.4301625 + 0.119727*x + -0.02510665*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.175636*x + -0.02611235*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.018007*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.260688 + 0.122812*x + -0.0242179*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.30177 + 0.136514*x + -0.02433725*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.5430035 + 0.2320075*x + -0.03007555*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.6874735 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.3817265 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.146904*x + 0.0156055*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.293806 + -0.1413035*x + 0.0129526*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0120754*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.01262435*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.4265315 + 0.20007*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.0139466*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.012562*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.763295 + -0.123621*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.8738605 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.5635785 + -0.161457*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.295472*x + 0.01852775*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.612026 + -0.1203225*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.429258 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.186897*x + 0.0164355*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.238211*x + 0.018642*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.263653 + -0.1259265*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.08729235*x + 0.008652835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0.3241825 + -0.22296*x + 0.0110055*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.096995 + -0.36493*x + 0.01945115*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.879221 + 0*x + -0.00371706*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.3347215*x + 0.0105112*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.4786015 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.8939105 + -0.3195775*x + 0.0142767*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.3411815 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.912655 + -0.2791495*x + 0.0117167*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0.315848 + -0.21914*x + 0.0152463*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-1.11077 + -0.5381785*x + 0.0308385*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0.3988775 + -0.092795*x + 0.008477605*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.576031 + -0.531151*x + 0.0279142*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + -0.09979145*x + 0.01171585*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.996043 + -0.424857*x + 0.02612505*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0.7338375 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6194855 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4650165 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4901915 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.194946*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.196648*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2148725*x + 0.024929*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(1.19861 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.72591 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.506579 + -0.1493705*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.393202*x + 0.0170683*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.595457 + -0.3009115*x + 0.0143427*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4860145 + -0.3033635*x + 0.01320535*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.61224 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.4030275*x + 0.01865645*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0.876623 + -0.1448815*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.3036135 + -0.1663285*x + 0.01991955*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.4556465 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.893063 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.381272 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.523034 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.391182 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7323775 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.114691*x + -0.0109624*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.008364605*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4003325 + 0*x + -0.005380905*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.645894 + 0.0881441*x + -0.01144835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.502073 + 0*x + -0.007972465*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.08179965*x + -0.0100508*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1055805*x + -0.00927189*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.006155595*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.487423 + 0*x + -0.00951779*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4596655 + 0.09921275*x + -0.01185675*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.382166 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){(0.2388755 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1254415*x + -0.01086365*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.006579255*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.583547 + 0.1826845*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7541255 + 0.1476565*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5309975 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.4414055 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.6630925 + 0.141492*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6900255 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5400705 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(-0.2324465 + -0.1726115*x + 0.007361315*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.353806 + -0.137524*x + 0.019214*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0.3698535 + -0.264*x + 0.02425245*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.3785625 + 0.205022*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2706855*x + 0.01345105*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.195684*x + 0.01244695*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.168131*x + 0.01204785*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.390939*x + 0.0255996*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(-0.295401 + -0.268837*x + 0.0175611*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1754695*x + 0.01133415*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.110268*x + 0.005343465*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1158865*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0.292162 + 0.1360315*x + -0.0110057*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.430456 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8523075 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.388355 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.3302575 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.0912052*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7341665 + 0.1477015*x + -0.01509915*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.753422 + 0.0823276*x + -0.01023145*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.399607 + 0*x + -0.0087317*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1440835*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.256078*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.528627 + 0.110531*x + -0.01789425*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1247235*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.264198 + -0.169129*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2675145*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.4428745 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1392025*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.969191 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.03897 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8844005 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(1.04319 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.7314515 + -0.236585*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.9770775 + -0.273102*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.6419825 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3271175 + 0*x + 0.032273*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.1791535*x + 0.0332422*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.1807145*x + 0.03567035*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.6882725 + -0.284502*x + 0.0367779*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.455788 + -0.2655285*x + 0.0350543*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.41844 + -0.170852*x + 0.03119915*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.164082*x + 0.02606895*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2727295*x + 0.04042885*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0.584382 + -0.188847*x + 0.03140805*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1173195*x + 0.0313835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1929605*x + 0.0314701*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.02833215*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1944185*x + 0.0317326*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2468965*x + 0.03795895*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.4942965 + -0.281549*x + 0.0391019*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.513104 + -0.2338995*x + 0.0350707*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.270964*x + 0.04065335*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.203775*x + 0.03586595*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0.278338 + -0.126748*x + 0.02956835*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.422215 + -0.1987455*x + 0.03202525*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.458044 + -0.159117*x + 0.0295487*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.056555 + -0.1409495*x + 0.0246246*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.624395 + 0.4319585*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.861475 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.616545 + 0.453692*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.34237 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1920865*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4182275 + -0.3117205*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.443508 + -0.351106*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.4424085*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.699695 + -0.1653025*x + -0.01623275*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2792245*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.237696*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(1.04209 + 0*x + -0.03017135*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.184705*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1693085*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2515595 + -0.1695895*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.6727815 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1813055*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.11825*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.212907*x + 0.0108101*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + -0.2118165*x + 0.01043015*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0.2563115 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.463193 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(1.023795 + 0*x + -0.01560605*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5338895 + 0.1060645*x + -0.01198095*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.4143605 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.466248 + -0.137611*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.822188 + 0*x + -0.0109008*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1262035*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.61744 + 0*x + -0.0157012*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.249551 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.331146 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.259762 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.4769895 + 0.190869*x + -0.021259*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(0.725389 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.7658975 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.5951095 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.6706775 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.74666 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.561176 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.4931045 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.271764 + -0.2500135*x + 0.0298893*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0.45494 + -0.151925*x + 0.02053925*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.243432*x + 0.01447565*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(1.07507 + 0.127439*x + -0.01300595*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1648525*x + 0.008554625*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0.850658 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.558802 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0.5918315 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0.335255 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0.2388595 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1975625*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.2627425 + -0.1797775*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.6819755 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.136475*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.480287 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.598073 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.784434 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.759163 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.8480005 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1219325*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.124214*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.2006875*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.4459 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.2660885 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.688117 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.637682 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.846348 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.6611805 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.7063845 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.507344 + -0.1736715*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(1.20301 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(0.3787 + -0.152318*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.3530725 + 0.2019165*x + -0.0222162*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.820268 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.522843 + 0.1502555*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.462557 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0.423135 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.5255085 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.37154 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.3708535 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.3131995 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.3884295 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(-0.976253 + 0*x + 0.01671*x^2)*(0.2287037)+(-0.07758351)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
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
                                  


# richnessInsetPlot <- function() {
#   print(richnessPlot)
#   
#   theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
#   
#   print(richnessInset, vp=viewport(width = 0.13, height = 0.25, x = 0.30, y = 0.99, just = c("left","top")))
#   
#   theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
# }

print(richnessPlot) #export at 1200x1000






###evenness
chainsEvennessIntercept <- chainsCommunity[,587:876]%>%
  gather(key=parameter, value=value, B.1.3.1:B.290.3.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), intercept=median(value))%>%
  mutate(lower=intercept-sd, upper=intercept+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, intercept=ifelse(diff==-2, 0, intercept))%>%
  select(parameter, intercept)
names(chainsEvennessIntercept)[1] <- 'parameter1'

chainsEvennessSlope <- chainsCommunity[,1747:2036]%>%
  gather(key=parameter, value=value, B.1.3.2:B.290.3.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), slope=median(value))%>%
  mutate(lower=slope-sd, upper=slope+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, slope=ifelse(diff==-2, 0, slope))%>%
  select(parameter, slope)
names(chainsEvennessSlope)[1] <- 'parameter2'

chainsEvennessQuad <- chainsCommunity[,2907:3196]%>%
  gather(key=parameter, value=value, B.1.3.3:B.290.3.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), quad=median(value))%>%
  mutate(lower=quad-sd, upper=quad+sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, quad=ifelse(diff==-2, 0, quad))%>%
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
  mutate(yr10=intercept + 10*slope + (10^2)*quad,
         yr20=intercept + 20*slope + (20^2)*quad,
         final_year_estimate=intercept + alt_length*slope + (alt_length^2)*quad)%>%
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




# #inset
# evennessInset <- ggplot(data=evenness, aes(x=final_year_estimate)) +
#   geom_density() +
#   xlab('Evenness Change') +
#   ylab('Density')


#main figure
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.35,0.55))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.2869265*x + -0.0165641*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2438115*x + -0.01359255*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.648768 + 0.171627*x + -0.01404535*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.6159045 + 0*x + -0.01400295*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.2999065 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.1238505*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6039485 + -0.3012245*x + 0.01349885*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.008677635*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1458705*x + -0.02350035*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1295205*x + -0.0225283*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0205062*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.2670335 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.348193 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.2938475 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.260833 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2690125 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.2728875 + -0.2106535*x + 0.01126395*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(-0.8076745 + -0.3056815*x + 0.01889215*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4594165 + -0.115953*x + 0.00640046*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.124157*x + 0.02155395*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.145456*x + 0.0244246*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0170701*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0.0202809*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4033845 + -0.250867*x + 0.0298407*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.128357*x + 0.02212545*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.562218 + 0.1116835*x + -0.0171216*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.659341 + 0*x + -0.01410315*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.891012 + 0*x + -0.0111684*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.667963 + 0.1111865*x + -0.0173492*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.8423015 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.012869*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.01125415*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0127721*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.01223115*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.01229045*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0131843*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0133048*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.01137005*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0136476*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.104102*x + -0.0139801*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1297325*x + -0.0178692*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3582955 + 0*x + -0.015439*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.164575*x + -0.020562*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.183452*x + -0.0232075*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.3224445 + 0*x + -0.014849*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2131235*x + -0.0164914*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.2085045*x + -0.01640775*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(-0.361842 + 0.124871*x + -0.00782899*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1086695*x + -0.00425099*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0.49791 + 0.486815*x + -0.0296503*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(-0.321141 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.342991*x + -0.01215835*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.3129785*x + -0.0142065*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){(0.2783285 + 0*x + -0.004478085*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.1693845*x + -0.0084182*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(1.40602 + 0.9065815*x + -0.06015*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(-0.3974515 + 0.0966749*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.630119*x + -0.0334193*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.07877575*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.4648995*x + -0.02317545*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4761505 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.431552 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(-0.2600835 + 0.156199*x + -0.0174344*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.197942*x + -0.0206162*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.319995 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.2750345 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3788885 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.01278025*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.0127801*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.266442 + 0.23086*x + -0.01991335*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.131717*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.474218 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
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
  stat_function(fun=function(x){(0.363212 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.271579 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.10059*x + 0.004992565*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.117303*x + 0.004452295*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.1067495*x + -0.00984472*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0103406*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6513185 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6305825 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.00528782*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){(0.296672 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0.2863585 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.08024575*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.1118045*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2850775 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.0778151*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.132035*x + 0.00599123*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.142919*x + 0.01092495*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1801455*x + -0.01036115*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.342502 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1753035*x + -0.0123585*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.184255*x + -0.0109484*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1356835*x + -0.007450275*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.0660244*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1163515*x + -0.00722828*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){(0.2536535 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0.361549 + 0.1760245*x + -0.0148833*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.331955 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4620645 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3074985 + 0.167376*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0.3367855 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.37288 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.5227815 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6155435 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.479389 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.333139 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3864445 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(-0.2886325 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.3881525 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.179416*x + -0.0141826*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.121826*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(-0.340264 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + -0.01254035*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + -0.01251065*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0108809*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.12949*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.4654085 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0 + -0.138861*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.294016 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-1.359805 + 0.1985735*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.33328 + 0.1388595*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-1.30181 + 0.126135*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2072825*x + 0.01374975*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.1978555*x + 0.01335345*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3358715 + 0.295233*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.523093*x + -0.0257636*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.383265 + 0.574114*x + -0.0297145*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0.403603 + 0.862082*x + -0.04261665*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.2804135 + 0.3199705*x + -0.0137251*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.3661685*x + -0.01533315*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(-0.306236 + 0.3524815*x + -0.01439955*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.268978 + 0.334261*x + -0.0128164*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.364569*x + -0.0135943*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.412672*x + -0.0172374*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.4349475*x + -0.01876915*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.1548045*x + -0.01214775*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2257315*x + -0.01506585*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(0 + 0.2377275*x + -0.01498675*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.2512495*x + -0.01643755*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0.203584*x + -0.01285325*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0.157275*x + -0.009493085*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0.01581605*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.4702325*x + -0.02763595*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0.0164506*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.266361*x + 0.02586255*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3298555*x + 0.0300639*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.208304*x + 0.0210999*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.2766255 + -0.3190975*x + 0.03012865*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){(-0.4077825 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.263883 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.256138 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.6836215 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6673645 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(-0.60223 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.484957 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){(-0.306087 + 0.132585*x + -0.0146953*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.177271*x + -0.01748665*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.114613*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.134768*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.538797 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4039595 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.132582*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(-0.3059935 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.2634545 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.104646*x + -0.00855397*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.217861*x + -0.0109906*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4987835 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0.08720315*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0.5550905 + 0.21163*x + -0.01412815*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(-0.285684 + 0.109451*x + -0.007786945*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.0908186*x + -0.007861805*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3328435 + 0.172183*x + -0.01827825*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4080185 + 0.1163005*x + -0.013911*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2041825*x + -0.0215562*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.295979 + 0.208523*x + -0.0203584*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.4502975 + 0.2114845*x + -0.02102815*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.509345 + 0.1167755*x + -0.0149689*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){(-0.398435 + 0.150529*x + -0.0192539*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.7058875 + 0*x + -0.01383875*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(-0.4895105 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.227075*x + 0.01910705*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + -0.3098325*x + 0.02541585*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0.263756 + -0.2699515*x + 0.0191664*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + -0.2845285*x + 0.0234434*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.3838405 + 0*x + 0.0134365*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5308165 + 0*x + -0.0134331*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0.2316735*x + -0.02039595*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(-0.340094 + 0.203116*x + -0.017684*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.6127565 + 0.2073355*x + -0.01652495*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(-0.5366045 + 0.5073425*x + -0.03066115*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){(-0.305906 + 0.2481755*x + -0.0209015*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){(-0.594625 + 0.471253*x + -0.0282798*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0.268412 + 0.2118155*x + -0.02539825*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2638215 + 0.158283*x + -0.0240929*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.01767595*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){(0 + 0*x + -0.0162135*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){(0.2833325 + 0.170015*x + -0.02535545*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){(0.269429 + 0.168616*x + -0.02475015*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.41435 + 0.2472745*x + -0.02915045*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.0152823*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0.3224535 + 0.1840225*x + -0.0242973*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + -0.019153*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){(0 + 0*x + -0.01295535*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){(0 + 0*x + -0.0187921*x^2)*(0.1034254)+(0.019179)}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  
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



# evennessInsetPlot <- function() {
#   print(evennessPlot)
#   
#   theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
#   
#   print(evennessInset, vp=viewport(width = 0.13, height = 0.25, x = 0.80, y = 0.99, just = c("left","top")))
#   
#   theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
# }

print(evennessPlot) #export at 1200x1000









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
