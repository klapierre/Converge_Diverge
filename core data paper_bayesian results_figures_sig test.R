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

chainsCommunity <- read.csv('fullChains.csv')

chainsMeanIntercept <- chainsCommunity[,7:296]%>%
  gather(key=parameter, value=value, B.1.1.1:B.290.1.1)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsMeanSlope <- chainsCommunity[,1167:1456]%>%
  gather(key=parameter, value=value, B.1.1.2:B.290.1.2)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

chainsMeanQuad <- chainsCommunity[,2327:2616]%>%
  gather(key=parameter, value=value, B.1.1.3:B.290.1.3)%>%
  group_by(parameter)%>%
  summarise(sd=sd(value), median=median(value))%>%
  mutate(lower=median-sd, upper=median+sd)

#mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
chainsCommunity2 <- chainsCommunity%>%
  select(lp__, U_int.2.1, U_int.2.2, U_int.2.3, U_int.2.4, U_slope.2.1, U_slope.2.2, U_slope.2.3, U_slope.2.4, U_quad.2.1, U_quad.2.2, U_quad.2.3, U_quad.2.4, mu_int.2, mu_slope.2, mu_quad.2, U_int.3.1, U_int.3.2, U_int.3.3, U_int.3.4, U_slope.3.1, U_slope.3.2, U_slope.3.3, U_slope.3.4, U_quad.3.1, U_quad.3.2, U_quad.3.3, U_quad.3.4, mu_int.3, mu_slope.3, mu_quad.3, U_int.1.1, U_int.1.2, U_int.1.3, U_int.1.4, U_slope.1.1, U_slope.1.2, U_slope.1.3, U_slope.1.4, U_quad.1.1, U_quad.1.2, U_quad.1.3, U_quad.1.4, mu_int.1, mu_slope.1, mu_quad.1, U_int.4.1, U_int.4.2, U_int.4.3, U_int.4.4, U_slope.4.1, U_slope.4.2, U_slope.4.3, U_slope.4.4, U_quad.4.1, U_quad.4.2, U_quad.4.3, U_quad.4.4, mu_int.4, mu_slope.4, mu_quad.4)%>%
  gather(key=parameter, value=value, U_int.2.1:mu_quad.4)%>%
  group_by(parameter)%>%
  summarise(median=median(value))

###mean change
mean <- read.csv('mean change_experiment_coefs.csv')%>%
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
# mean2 <- mean%>%
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


#inset - density plot of mean change (all datapoints)
meanInset <- ggplot(data=mean, aes(x=final_year_estimate)) +
  geom_density() +
  xlab('Mean Change') +
  ylab('Density')


#main figure
meanPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  ylim(-1,2) +
  xlab('Standardized Year') +
  ylab('Mean Change')

meanPlot <- meanPlot + 
#below are the individual treatment lines
  stat_function(fun=function(x){0.230929386 + 0.021502518*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.262171061 + 0.02401176*x + -0.001525823*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.264315375 + 0.030791937*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.206609096 + 0.031686224*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.227275001 + 0.021500391*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.270091532 + 0.027955621*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.248493062 + 0.029207605*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.25755068 + 0.029005406*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.034425311*x + -0.001183539*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.252928172 + 0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.260585623 + 0.04300529*x + -0.001497189*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.219919105 + 0.051229869*x + -0.0015955*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.37634329 + 0.056351027*x + -0.002208445*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.433885054 + 0.065284875*x + -0.002671376*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.530320075 + 0.054212157*x + -0.003370371*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.505870741 + 0*x + -0.001773763*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.488987074 + 0*x + -0.001834057*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.484917572 + 0.019967268*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.387016459 + 0.039793413*x + -0.002103177*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.41234162 + 0.035825734*x + -0.001983006*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.389915554 + 0.035089413*x + -0.001925995*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.370785749 + 0.024994854*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.042219887*x + -0.001975358*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.413064501 + 0.045393826*x + -0.00241283*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.375913542 + 0.036219329*x + -0.001965967*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.032617514*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.02973926*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.031268726*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.030973806*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.351048582 + 0.037485093*x + -0.001483422*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.176098128 + 0.023011823*x + -0.001517118*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.197645645 + 0.016056216*x + -0.001328941*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.139780634 + 0.018570418*x + -0.00128736*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.178500188 + 0.034384225*x + -0.001931856*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.21031503 + 0.04433613*x + -0.002585511*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.168717393 + 0.023587201*x + -0.001663977*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.171530146 + 0.012670687*x + -0.001056369*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.056500464 + 0.03138637*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.036331593 + 0.021532631*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.043494051 + 0.016490855*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.036586787 + 0.017710923*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.058007813 + 0.025308658*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.031204735 + 0.028832469*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.205895487 + 0*x + -0.002044338*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.228201187 + 0*x + -0.001977204*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.029693835*x + -0.00294958*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.02500285*x + -0.002734868*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.019490479*x + -0.002490953*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.21825915 + 0*x + -0.002228367*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.233517739 + 0.029103656*x + -0.002951409*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.209729359 + 0.015962118*x + -0.002418104*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.229964751 + 0.015391477*x + -0.002460041*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.250137875 + 0*x + -0.00226389*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.01690313*x + -0.002331721*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.256035165 + 0*x + -0.001919981*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.280592275 + 0.025900795*x + -0.002667701*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.300746635 + 0.031266174*x + -0.002934133*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.240956573 + 0.014305837*x + -0.002156138*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.246521514 + 0.031779115*x + -0.00242456*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.209583558 + 0*x + -0.001562067*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.273619681 + 0.029979143*x + -0.002424654*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.041813192*x + -0.00289906*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.200072034 + 0.017587664*x + -0.001831191*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.199834363 + 0.01524691*x + -0.001158153*x^2}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){0.193679582 + 0.021933711*x + -0.00151169*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){0.092756796 + 0*x + -0*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.414275994 + 0.054686563*x + -0.001778748*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.162265736 + 0*x + -0*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.238162534 + 0.090281176*x + -0.002742005*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.23110513 + -0*x + 0*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.431619097 + 0.034674722*x + -0.001197686*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.332504534 + 0.022716138*x + -0.001433731*x^2}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){0.526499813 + 0.051668293*x + -0.002729824*x^2}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){0.212155578 + 0*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.042198876*x + -0.001606892*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.274940653 + 0*x + -0.001127782*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.488647665 + 0.092501368*x + -0.004392884*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.224272382 + 0*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.354865101 + 0.079669083*x + -0.00290311*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.185333106 + 0.018662883*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.216815344 + 0.023667588*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.18501105 + 0.043916505*x + -0.001759745*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.134870692 + 0.029399766*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){0.084538683 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.103025822 + 0.02063256*x + -0.001783486*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.033612517*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.243143675 + 0.030087941*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.239141204 + 0.05581129*x + -0.001529445*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.170299003 + 0.038571371*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.231644356 + 0.045839225*x + -0.00168189*x^2}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.037862866*x + -0.001476611*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.033613708*x + -0^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.183128566 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){0.221359678 + 0.031098681*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.206548019 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.177079861 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.168546157 + 0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.160470783 + 0.022087083*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.111038079 + 0.019262931*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.141435145 + 0.02610461*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.126222151 + 0.020588836*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.122275143 + 0.019326985*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.148797591 + 0.028556604*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.213510916 + 0.032284485*x + -0.001890864*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.245896798 + 0.023463602*x + -0.001633979*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.239492862 + 0.029539103*x + -0.001843644*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.24326889 + 0.021690511*x + -0.001472961*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.249007619 + 0.025506519*x + -0.00170458*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.236549194 + 0.035860355*x + -0.002274106*x^2}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){0.225421778 + 0.032975041*x + -0.00194032*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.251787793 + 0.033789962*x + -0.002161701*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.225248501 + 0.033361235*x + -0.001924413*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.2261655 + 0.027389684*x + -0.001734446*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.243098165 + 0.030633377*x + -0.001811821*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.203956434 + 0.022703548*x + -0.001422709*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.213013031 + 0.019043464*x + -0.001372079*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.206904101 + 0.02252959*x + -0.001641347*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.387045892 + 0.03220648*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.351790432 + 0.032924428*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.351972216 + 0.026820346*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.02847979*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.357430231 + 0.02758678*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.360967481 + 0.02104657*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.377678382 + 0.027230103*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.034229322*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + 0*x^2}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){0.646820613 + 0.024231738*x + -0.000792657*x^2}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){0.171654001 + 0*x + -0*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.244530742 + 0.018780528*x + -0.001734863*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.196132001 + 0.023560746*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.160551169 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.402715514 + 0.018596363*x + -0.000976014*x^2}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.021489928*x + -0.000880253*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.356586047 + 0.033082393*x + -0.001035054*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.373853782 + 0.03420627*x + -0.001352192*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.378325215 + 0.041800687*x + -0.001948674*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.021443568*x + -0.000586184*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.386339088 + 0.013074915*x + -0.000804871*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.452625602 + 0.022434998*x + -0.001418926*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.445655219 + 0.033188043*x + -0.002082931*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.361731279 + 0.021366924*x + -0.001116561*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.010623636*x + -0*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0.193902622 + 0.029048534*x + -0.001235522*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.18373831 + 0.02216798*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.20742189 + 0.020333642*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.176285781 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.149509414 + 0.02430702*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.131342203 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.134903017 + 0.025878508*x + -0.000777135*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0.167058713 + 0*x + 0*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0.186414875 + 0*x + -0*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.216915636 + 0.036364619*x + -0.001814552*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.249747172 + 0.044406393*x + -0.002199385*x^2}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){0.175268831 + 0*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.142193073 + 0.012028609*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.161395778 + 0.028078369*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.143497116 + 0.03288198*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.20518477 + 0*x + 0.00155233*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.241345234 + 0.032013213*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.031379055*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + 0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.031664277*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.359924331 + 0.045986217*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.266884163 + 0.056284421*x + -0.001675281*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.280226497 + 0.053399277*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.279816399 + 0.036711939*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.296853353 + 0.038498471*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.23399189 + 0.032650008*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.167682834 + 0.02130891*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.177982399 + 0.025512814*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.181873519 + 0.025245795*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.246126728 + 0.029025311*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.250579107 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.165482207 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.251976637 + 0.01974712*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){0.215662461 + 0.016847421*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.228364766 + 0.014081113*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.243304873 + 0.043823614*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){0.253866777 + 0.042833119*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.217004273 + 0.039462425*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.179172881 + 0*x + 0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.176647902 + 0.033122629*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.17276163 + 0*x + 0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.203789707 + 0.021451819*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.1528175 + 0*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.197601496 + 0.017488478*x + -0*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.026167983*x + -0.001276188*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.02048846*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.023714714*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.256341823 + 0.031811269*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.029523195*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){0.252979381 + 0.038472101*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.154122394 + 0.027309213*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.209291785 + 0.043498156*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.206741457 + 0.028264916*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.208998907 + 0.025900114*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.026762927*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.226807314 + 0.018823826*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.037523202*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.237639555 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.259423382 + 0.0258683*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.359927819 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.432920844 + 0.044692722*x + -0.002807794*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.381273052 + 0.034273045*x + -0.001820549*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.019653038*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.02497614*x + -0.001950477*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.356919927 + 0.016319271*x + -0.001616396*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.018223354*x + -0.001606864*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.235814574 + 0*x + -0.001627748*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.270821389 + 0.02280435*x + -0.001915609*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.208666814 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.19018699 + 0*x + 0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.05430258*x + -0.001866135*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.047485995*x + -0.001655464*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.047808901*x + -0.001617279*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.053605219*x + -0.001918748*x^2}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){0.171317739 + 0.016972382*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.118306868 + 0.033496233*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.18628311 + 0.065670559*x + -0.001807713*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.124018121 + 0.030336841*x + -0*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.093148945 + 0.01411747*x + 0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.122500564 + 0.033656155*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.124303939 + 0.022475659*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.145399762 + 0.031374036*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){0.143675753 + 0.020362734*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.147866982 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.141501496 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.158699733 + 0*x + -0*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.166254255 + 0.016904781*x + -0.001413259*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.235672941 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){0.066290575 + 0*x + -0*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){0.199322613 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.224405849 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.221920595 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.186889707 + 0.018288768*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.189141799 + 0*x + -0*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.235925243 + 0.018755179*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.215343212 + 0.02148559*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.245393044 + 0.019511065*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.223195887 + 0.031605923*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.271063994 + 0*x + -0.001837775*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.00182886*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.197506308 + 0.023498734*x + -0.002490783*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.210232432 + 0.015816742*x + -0.002536284*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.077388133 + 0*x + -0.001474587*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.087784757 + 0*x + -0.001714932*x^2}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){0.158780715 + 0*x + 0*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.000840271*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.03341472 + 0.037417212*x + -0.002066182*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.22262187 + 0*x + -0.001475355*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.130798639 + 0.039854914*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.129730224 + 0.041991063*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.133538577 + 0.0286637*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.133221285 + 0.036650607*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.206606629 + 0.014993791*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.169838972 + 0.025115731*x + -0.001476843*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.185743373 + 0.037495216*x + -0.001913474*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.165812088 + 0.024629841*x + -0.001475891*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.037698096*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.220553263 + 0.023792208*x + 0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.046494139*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.15727932 + 0*x + 0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.215567698 + 0.02113793*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0.018634047*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.234500067 + 0.031099532*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.25887029 + 0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + -0*x + -0*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0*x + -0.001850075*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.440077859 + 0*x + -0.002206003*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.401007412 + 0*x + -0.002004315*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.39534907 + 0*x + -0.001860887*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0 + 0.028710656*x + -0.002300434*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0.00191133*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0 + 0.033906501*x + -0.002448744*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.162807428 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.255856784 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.208078846 + 0.019881012*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0 + 0*x + -0*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.322484612 + 0.007737378*x + -0.001270659*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.247561942 + 0.024541204*x + -0.001652742*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.241731854 + 0.031946777*x + -0.00197068*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.273621638 + 0.032687607*x + -0.002125387*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.028155183*x + -0.001922593*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.029392281*x + -0.002033441*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0.420676102 + 0.031447617*x + -0.002369625*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0 + 0.023773068*x + -0.001845754*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0.276706854 + 0*x + -0*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
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



meanInsetPlot <- function() {
  print(meanPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(meanInset, vp=viewport(width = 0.13, height = 0.25, x = 0.115, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

meanInsetPlot() #export at 1200x1000




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
