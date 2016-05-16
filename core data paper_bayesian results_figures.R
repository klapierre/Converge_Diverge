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
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon))

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
meanInset <- ggplot(data=rawData, aes(x=mean_change)) +
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
  stat_function(fun=function(x){0.230929386 + 0.021502518*x + -0.001364054*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.262171061 + 0.02401176*x + -0.001525823*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.264315375 + 0.030791937*x + -0.001286028*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.206609096 + 0.031686224*x + -0.001175413*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.227275001 + 0.021500391*x + -0.00082508*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.270091532 + 0.027955621*x + -0.001133344*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.248493062 + 0.029207605*x + -0.001166054*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.25755068 + 0.029005406*x + -0.001172893*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.285460706 + 0.034425311*x + -0.001183539*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.252928172 + 0.001981279*x + -0.000406796*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.260585623 + 0.04300529*x + -0.001497189*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.219919105 + 0.051229869*x + -0.0015955*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.37634329 + 0.056351027*x + -0.002208445*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.433885054 + 0.065284875*x + -0.002671376*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.530320075 + 0.054212157*x + -0.003370371*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.505870741 + 0.009981243*x + -0.001773763*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.488987074 + 0.013030401*x + -0.001834057*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.484917572 + 0.019967268*x + -0.001608783*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.387016459 + 0.039793413*x + -0.002103177*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.41234162 + 0.035825734*x + -0.001983006*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.389915554 + 0.035089413*x + -0.001925995*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.370785749 + 0.024994854*x + -0.001472578*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.310387372 + 0.042219887*x + -0.001975358*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.413064501 + 0.045393826*x + -0.00241283*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.375913542 + 0.036219329*x + -0.001965967*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.338939008 + 0.032617514*x + -0.001255627*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.287855876 + 0.02973926*x + -0.001099468*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.332686233 + 0.031268726*x + -0.001218994*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.350133114 + 0.030973806*x + -0.00123076*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.351048582 + 0.037485093*x + -0.001483422*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.176098128 + 0.023011823*x + -0.001517118*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.197645645 + 0.016056216*x + -0.001328941*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.139780634 + 0.018570418*x + -0.00128736*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.178500188 + 0.034384225*x + -0.001931856*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.21031503 + 0.04433613*x + -0.002585511*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.168717393 + 0.023587201*x + -0.001663977*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.171530146 + 0.012670687*x + -0.001056369*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.056500464 + 0.03138637*x + -0.001115602*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.036331593 + 0.021532631*x + -0.000853358*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.043494051 + 0.016490855*x + -0.000756066*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.036586787 + 0.017710923*x + -0.000888481*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.058007813 + 0.025308658*x + -0.001014379*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.031204735 + 0.028832469*x + -0.000926348*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.205895487 + 0.011412366*x + -0.002044338*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.228201187 + 0.008703459*x + -0.001977204*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.324029763 + 0.029693835*x + -0.00294958*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.297351986 + 0.02500285*x + -0.002734868*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.285805899 + 0.019490479*x + -0.002490953*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.21825915 + 0.009963729*x + -0.002228367*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.233517739 + 0.029103656*x + -0.002951409*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.209729359 + 0.015962118*x + -0.002418104*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.229964751 + 0.015391477*x + -0.002460041*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.250137875 + 0.0079407*x + -0.00226389*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.283103219 + 0.01690313*x + -0.002331721*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.256035165 + 0.00629114*x + -0.001919981*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.280592275 + 0.025900795*x + -0.002667701*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.300746635 + 0.031266174*x + -0.002934133*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.240956573 + 0.014305837*x + -0.002156138*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.246521514 + 0.031779115*x + -0.00242456*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.209583558 + 0.009936593*x + -0.001562067*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.273619681 + 0.029979143*x + -0.002424654*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.298308124 + 0.041813192*x + -0.00289906*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.200072034 + 0.017587664*x + -0.001831191*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.199834363 + 0.01524691*x + -0.001158153*x^2}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){0.193679582 + 0.021933711*x + -0.00151169*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){0.092756796 + 0.00947595*x + -0.000472847*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){0.300178112 + 0.004406443*x + -0.000272609*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.414275994 + 0.054686563*x + -0.001778748*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.162265736 + 0.007515903*x + -0.000300267*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.238162534 + 0.090281176*x + -0.002742005*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.23110513 + -0.001098398*x + 0.000290017*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.431619097 + 0.034674722*x + -0.001197686*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.332504534 + 0.022716138*x + -0.001433731*x^2}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){0.526499813 + 0.051668293*x + -0.002729824*x^2}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){0.212155578 + 0.009646862*x + -0.000636513*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.297137079 + 0.042198876*x + -0.001606892*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.274940653 + 0.00633864*x + -0.001127782*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.488647665 + 0.092501368*x + -0.004392884*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.224272382 + 0.010444898*x + -0.000534399*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.354865101 + 0.079669083*x + -0.00290311*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.185333106 + 0.018662883*x + -0.001116384*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.216815344 + 0.023667588*x + -0.001387654*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.18501105 + 0.043916505*x + -0.001759745*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.134870692 + 0.029399766*x + -0.001019035*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){0.084538683 + 0.010773588*x + -0.001416068*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.103025822 + 0.02063256*x + -0.001783486*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.289375134 + 0.033612517*x + -0.000861912*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.243143675 + 0.030087941*x + -0.000573828*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.239141204 + 0.05581129*x + -0.001529445*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.170299003 + 0.038571371*x + -0.001190205*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.231644356 + 0.045839225*x + -0.00168189*x^2}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){0.300760687 + 0.037862866*x + -0.001476611*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.330861924 + 0.033613708*x + -0.001375621*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.183128566 + 0.006690919*x + -0.000390798*x^2}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){0.221359678 + 0.031098681*x + -0.001162498*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.206548019 + 0.015010719*x + -0.00032087*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.177079861 + 0.012995269*x + -5.44e-05*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.168546157 + 0.010486112*x + 2.88e-05*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.160470783 + 0.022087083*x + -0.000720115*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.111038079 + 0.019262931*x + -0.000642474*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.141435145 + 0.02610461*x + -0.000973803*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.126222151 + 0.020588836*x + -0.000727862*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.122275143 + 0.019326985*x + -0.000683863*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.148797591 + 0.028556604*x + -0.001083467*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
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
  stat_function(fun=function(x){0.387045892 + 0.03220648*x + -0.001141764*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.351790432 + 0.032924428*x + -0.001095016*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.351972216 + 0.026820346*x + -0.000793766*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.337911425 + 0.02847979*x + -0.000839261*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.357430231 + 0.02758678*x + -0.00088651*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.360967481 + 0.02104657*x + -0.000604602*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.377678382 + 0.027230103*x + -0.000919029*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.344820816 + 0.034229322*x + -0.001125056*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.298235368 + -0.000747876*x + 0.000243963*x^2}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){0.646820613 + 0.024231738*x + -0.000792657*x^2}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){0.171654001 + 0.005861324*x + -0.001018593*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.244530742 + 0.018780528*x + -0.001734863*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.196132001 + 0.023560746*x + -0.000659286*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.160551169 + 0.013483499*x + -0.000241208*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.402715514 + 0.018596363*x + -0.000976014*x^2}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){0.345973614 + 0.021489928*x + -0.000880253*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.356586047 + 0.033082393*x + -0.001035054*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.373853782 + 0.03420627*x + -0.001352192*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.378325215 + 0.041800687*x + -0.001948674*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.346667233 + 0.021443568*x + -0.000586184*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.386339088 + 0.013074915*x + -0.000804871*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.452625602 + 0.022434998*x + -0.001418926*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.445655219 + 0.033188043*x + -0.002082931*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.361731279 + 0.021366924*x + -0.001116561*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.337818449 + 0.010623636*x + -0.000257857*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0.301573124 + 0.010723715*x + -0.000383004*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0.193902622 + 0.029048534*x + -0.001235522*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.18373831 + 0.02216798*x + -0.000930215*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.20742189 + 0.020333642*x + -0.001039395*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.176285781 + 0.013252862*x + -0.000687049*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.149509414 + 0.02430702*x + -0.001013057*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.131342203 + 0.015690021*x + -0.000666863*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.134903017 + 0.025878508*x + -0.000777135*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0.167058713 + 0.004461523*x + 0.000105696*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0.186414875 + 0.006639574*x + -0.000399876*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.216915636 + 0.036364619*x + -0.001814552*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.249747172 + 0.044406393*x + -0.002199385*x^2}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){0.175268831 + 0.005226502*x + -0.000395803*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.142193073 + 0.012028609*x + 0.000137172*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.161395778 + 0.028078369*x + -0.000343627*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.143497116 + 0.03288198*x + -0.000527793*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.20518477 + 0.00015646*x + 0.00155233*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.241345234 + 0.032013213*x + 0.000496587*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.31587413 + 0.031379055*x + 0.000320753*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.314182991 + 0.003520255*x + 0.000291995*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.32631611 + 0.031664277*x + -0.000511443*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.359924331 + 0.045986217*x + -0.000865623*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.266884163 + 0.056284421*x + -0.001675281*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.280226497 + 0.053399277*x + -0.001549594*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.279816399 + 0.036711939*x + -0.001209437*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.296853353 + 0.038498471*x + -0.00123178*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.23399189 + 0.032650008*x + -0.001264682*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.167682834 + 0.02130891*x + -0.000699754*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.177982399 + 0.025512814*x + -0.000869724*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.181873519 + 0.025245795*x + -0.000866869*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.246126728 + 0.029025311*x + -0.001080994*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.250579107 + 0.016959069*x + -0.000719076*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.165482207 + 0.016323005*x + -0.000525092*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.251976637 + 0.01974712*x + -0.000780308*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){0.215662461 + 0.016847421*x + -0.000734618*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.228364766 + 0.014081113*x + -0.00065569*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.243304873 + 0.043823614*x + -0.000778101*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){0.253866777 + 0.042833119*x + -0.000772412*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.217004273 + 0.039462425*x + -0.000530521*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.179172881 + 0.007500438*x + 0.000490577*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.176647902 + 0.033122629*x + -0.000405396*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.17276163 + 0.009962308*x + 6.57e-05*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.203789707 + 0.021451819*x + -0.000184531*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.1528175 + 0.003900707*x + -0.000262252*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.197601496 + 0.017488478*x + -0.000720524*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.283507447 + 0.026167983*x + -0.001276188*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.290673479 + 0.016324596*x + -0.000258838*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.334444182 + 0.02048846*x + -0.000486821*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.338435254 + 0.023714714*x + -0.00062812*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.256341823 + 0.031811269*x + -0.000669143*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){0.299758658 + 0.029523195*x + -0.00065442*x^2}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){0.252979381 + 0.038472101*x + -0.000924434*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.154122394 + 0.027309213*x + -0.000447056*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.209291785 + 0.043498156*x + -0.001206248*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.206741457 + 0.028264916*x + -0.000591282*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.208998907 + 0.025900114*x + -0.001134547*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.336826168 + 0.026762927*x + -0.001416737*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.226807314 + 0.018823826*x + -0.000883024*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.294033029 + 0.037523202*x + -0.001579189*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.237639555 + 0.014729341*x + -0.000608712*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.259423382 + 0.0258683*x + -0.001655734*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.359927819 + 0.008296517*x + -0.001248443*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.432920844 + 0.044692722*x + -0.002807794*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.381273052 + 0.034273045*x + -0.001820549*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.315438784 + 0.019653038*x + -0.001562072*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.345452932 + 0.02497614*x + -0.001950477*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.356919927 + 0.016319271*x + -0.001616396*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.303416411 + 0.018223354*x + -0.001606864*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.235814574 + 0.01695872*x + -0.001627748*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.270821389 + 0.02280435*x + -0.001915609*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.208666814 + 0.013038006*x + -0.001426277*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.19018699 + 0.011453367*x + 2.42e-05*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.305151708 + 0.05430258*x + -0.001866135*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){0.316049381 + 0.047485995*x + -0.001655464*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0.305648989 + 0.047808901*x + -0.001617279*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0.340442699 + 0.053605219*x + -0.001918748*x^2}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){0.171317739 + 0.016972382*x + -0.000176419*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.118306868 + 0.033496233*x + -0.000576232*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.18628311 + 0.065670559*x + -0.001807713*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.124018121 + 0.030336841*x + -0.000774999*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.093148945 + 0.01411747*x + 3.03e-05*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.122500564 + 0.033656155*x + -0.00076154*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.124303939 + 0.022475659*x + -0.000395175*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.145399762 + 0.031374036*x + -0.000778897*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){0.143675753 + 0.020362734*x + -0.001355726*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.147866982 + 0.00952345*x + -0.0012212*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.141501496 + 0.00538861*x + -0.001018939*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.158699733 + 0.011054617*x + -0.001206203*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.166254255 + 0.016904781*x + -0.001413259*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.235672941 + 0.010639024*x + -0.00131296*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){0.066290575 + 0.012227899*x + -0.001065694*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){0.199322613 + 0.013967313*x + -0.001127102*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.224405849 + 0.016524498*x + -0.001285604*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.221920595 + 0.010422075*x + -0.000595084*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.186889707 + 0.018288768*x + -0.000946616*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.189141799 + 0.015973014*x + -0.000896071*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.235925243 + 0.018755179*x + -0.001116732*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.215343212 + 0.02148559*x + -0.001145399*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.245393044 + 0.019511065*x + -0.001142569*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.223195887 + 0.031605923*x + -0.00144141*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.271063994 + 0.009853774*x + -0.001837775*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.304939735 + 0.009329587*x + -0.00182886*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.197506308 + 0.023498734*x + -0.002490783*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.210232432 + 0.015816742*x + -0.002536284*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.077388133 + 0.007560247*x + -0.001474587*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.087784757 + 0.010177828*x + -0.001714932*x^2}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){0.158780715 + 0.000138247*x + 0.000415305*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.282946615 + 0.008871207*x + -0.000840271*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.03341472 + 0.037417212*x + -0.002066182*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.22262187 + 0.007107566*x + -0.001475355*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.130798639 + 0.039854914*x + -0.000761886*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.129730224 + 0.041991063*x + -0.000812829*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.133538577 + 0.0286637*x + -0.000314958*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.133221285 + 0.036650607*x + -0.000616581*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.206606629 + 0.014993791*x + -0.001190588*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.169838972 + 0.025115731*x + -0.001476843*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.185743373 + 0.037495216*x + -0.001913474*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.165812088 + 0.024629841*x + -0.001475891*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.322452074 + 0.037698096*x + -0.000725417*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.220553263 + 0.023792208*x + 7.84e-05*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.28508591 + 0.046494139*x + -0.000930989*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.15727932 + 0.011402923*x + 1.12e-05*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.215567698 + 0.02113793*x + -0.000511491*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.282896682 + 0.018634047*x + -0.000471229*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.234500067 + 0.031099532*x + -0.000939184*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.25887029 + 0.015193302*x + -0.001467677*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.314332152 + -0.001159635*x + -0.00096307*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.300114833 + 0.007972803*x + -0.001850075*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.440077859 + 0.008953618*x + -0.002206003*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.401007412 + 0.006314294*x + -0.002004315*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.39534907 + 0.002591441*x + -0.001860887*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.348203078 + 0.028710656*x + -0.002300434*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.327948394 + 0.006889469*x + -0.00191133*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.322579102 + 0.033906501*x + -0.002448744*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.162807428 + 0.01033929*x + -0.000990086*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.255856784 + 0.015100504*x + -0.00142521*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.208078846 + 0.019881012*x + -0.00144848*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.256094115 + 0.020393187*x + -0.001511413*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.322484612 + 0.007737378*x + -0.001270659*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.247561942 + 0.024541204*x + -0.001652742*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.241731854 + 0.031946777*x + -0.00197068*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.273621638 + 0.032687607*x + -0.002125387*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.298236346 + 0.028155183*x + -0.001922593*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.28492165 + 0.029392281*x + -0.002033441*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0.420676102 + 0.031447617*x + -0.002369625*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.297726884 + 0.023773068*x + -0.001845754*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0.276706854 + 0.011850926*x + -0.001318645*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#last five are the main plot_mani effect lines
#estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){((-0.5441655 + 0.18185) + (0.132175+0.03107995)*x + (-0.00629647-0.00158602)*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.5441655 + 0.132175*x + -0.00629647*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.5441655 + 0.731963) + (0.132175 + 0.212783)*x + (-0.00629647-0.007026515)*x^2)*0.1701297 + 0.3140121}, size=3, xlim=c(0,22), colour='#EC1804')





meanInsetPlot <- function() {
  print(meanPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(meanInset, vp=viewport(width = 0.13, height = 0.25, x = 0.115, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

meanInsetPlot() #export at 1200x1000




###dispersion
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
dispersionInset <- ggplot(data=rawData, aes(x=dispersion_change)) +
  geom_density() +
  xlab('Dispersion Change') +
  ylab('Density')


#main figure
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,24), ylim=c(-0.3,0.3))  +
  scale_x_continuous(limits=c(0,24), breaks=seq(1,24,2)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){0.010688636 + -0.011751124*x + -9.85e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.017571786 + -0.015890051*x + 0.000160694*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.007977054 + -0.011468854*x + 0.00022684*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.005325249 + -0.012895299*x + 0.001485225*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.02191456 + -0.014966281*x + 0.001482447*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.001085616 + -0.019016194*x + 0.001742704*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.064169496 + -0.017841245*x + 0.00159046*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.048682274 + -0.016542292*x + 0.001607859*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00021065 + 0.015905597*x + 0.000359621*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.020223806 + -0.008730244*x + 0.001327016*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.002154098 + -0.01386412*x + 0.001440809*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.014368443 + -0.007976371*x + 0.001283393*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.013761132 + -0.013431967*x + 0.001598169*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.006260501 + -0.019460086*x + 0.001825246*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.02933066 + -0.01089552*x + -0.000624543*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.037927391 + 0.004481976*x + -0.00136027*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.022398925 + -9.57e-05*x + -0.001206548*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027245251 + -0.028956762*x + 0.004383552*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.048314479 + -0.025357947*x + 0.004179582*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.060453159 + -0.02624632*x + 0.004076132*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.03776382 + -0.023438616*x + 0.004032436*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.092039053 + -0.023328118*x + 0.004064221*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026244749 + -0.034279431*x + 0.004698383*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016939834 + -0.030599171*x + 0.00448058*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.018613652 + -0.029626316*x + 0.004444485*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025015276 + 0.002265009*x + 0.000946726*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029885487 + 0.001969644*x + 0.00089501*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.062337365 + 0.003678238*x + 0.000790774*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.003083587 + 0.002239692*x + 0.001018246*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.034006512 + 0.00285998*x + 0.00091166*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.03492158 + -0.004670034*x + 0.000272997*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046096515 + 0.002290897*x + 6.24e-05*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00030435 + -0.010735304*x + 0.000733447*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00199633 + -0.004355856*x + 0.00038454*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.030383691 + -0.00487792*x + 0.000382396*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00631454 + -0.004371397*x + 0.000270136*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00725764 + -0.004933264*x + 0.000231877*x^2}, size=0.5, xlim=c(0,15), colour='#bdbdbd') +
  stat_function(fun=function(x){0.002166886 + 0.001600209*x + 0.00016932*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.005453655 + -0.002340476*x + 0.000287673*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.017902552 + -0.005397782*x + 0.000393133*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.01593903 + 0.003618988*x + 0.000201706*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.011773665 + -0.001154232*x + 0.000217684*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.000617104 + -0.005456408*x + 0.000351817*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.031832979 + -0.013380707*x + 0.002480388*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.02371303 + -0.013567029*x + 0.002453951*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.000918401 + -0.010890036*x + 0.002277486*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.03382093 + -0.014920641*x + 0.00265144*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.044805192 + -0.016957948*x + 0.002673399*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.005625164 + -0.000233747*x + 0.001999204*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.040278664 + -0.007524752*x + 0.002435196*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.019238216 + -0.012802342*x + 0.002587254*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.069755551 + -0.005880221*x + 0.002555515*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.034333123 + -0.008185051*x + 0.00248293*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.041189925 + -0.015925992*x + 0.002650321*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.045045358 + -0.014105102*x + 0.002588786*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.02865109 + -0.008989767*x + 0.002345792*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.068277392 + -0.007173305*x + 0.002489942*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.046407853 + -0.016131033*x + 0.002659608*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009962972 + -0.001362667*x + 0.002196522*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.033029366 + -0.006627119*x + 0.002537997*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00243274 + -0.009979001*x + 0.002596487*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.010416971 + -0.002150156*x + 0.002330419*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.000401332 + -0.002988121*x + 0.002316246*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.005178228 + 0.004919309*x + -8.02e-05*x^2}, size=0.5, xlim=c(0,13), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.004205849 + 0.004428155*x + -0.000201971*x^2}, size=0.5, xlim=c(0,13), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026248964 + -0.000937063*x + -8.28e-05*x^2}, size=0.5, xlim=c(0,13), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.080184064 + 0.004993915*x + 8.87e-06*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.210704821 + 0.003895824*x + 0.000188387*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.004255278 + -0.007903859*x + 0.000451697*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.096746283 + -0.015661624*x + 0.000118136*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.007695788 + -0.011028361*x + 0.000580176*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.066514953 + 0.017938236*x + -0.000965136*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){0.036551434 + 0.003794292*x + -0.000246995*x^2}, size=0.5, xlim=c(0,20), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.001066989 + 0.006134401*x + -0.000382569*x^2}, size=0.5, xlim=c(0,20), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.034844456 + -0.004074559*x + -0.000236525*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.139358515 + -0.026485897*x + 0.000712483*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00120967 + -0.001713376*x + -0.000158627*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.015120998 + -0.035370941*x + 0.000899866*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.013338361 + 0.006469804*x + -0.000125287*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.043426773 + 0.00716424*x + -0.000621176*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.013887114 + -0.003276651*x + 0.000377638*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027205729 + -0.003416771*x + 0.000452866*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.002383702 + 0.003316848*x + 0.001015259*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017149679 + 9.38e-05*x + 0.001157432*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.001995788 + -0.003164658*x + 0.001921788*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009126494 + 0.006648512*x + 0.001605906*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.027009088 + -0.008989626*x + 0.002532903*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.01107342 + -0.006186228*x + 0.002379028*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.025818593 + -0.002717489*x + 0.002399083*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.009368823 + -0.013208843*x + 0.001823152*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.00141804 + -0.010660158*x + 0.001535579*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.011750233 + -0.022274045*x + 0.000535022*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024704815 + -0.024483534*x + 0.000662006*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.012731337 + -0.023498986*x + 0.000510167*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.054257708 + -0.044595679*x + 0.001238737*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.108294996 + 0.009690793*x + -0.000489199*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.005699717 + 0.003328949*x + -9.69e-05*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.026781568 + 0.00467409*x + -9.74e-05*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.078678711 + -0.005380238*x + 0.000196923*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.027245674 + 0.004109562*x + -7.33e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.006536028 + 0.00544238*x + -0.000111388*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.027514801 + 0.005795377*x + -0.00018123*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.054222689 + 0.004442644*x + -4.64e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.016255474 + 0.003091226*x + 6.95e-06*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.030289238 + 0.001064484*x + 1.78e-05*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.010842507 + -0.002237602*x + 0.00011896*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024459709 + -0.001491865*x + 0.000121224*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.012896992 + -0.004184209*x + 0.000286503*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027250281 + 0.002256778*x + -0.000205675*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.006541276 + -0.004040368*x + 0.000249605*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.028448572 + -0.006448284*x + 0.000213771*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.000905445 + -0.002205767*x + 0.000206285*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.004035139 + -0.004452647*x + 0.000180665*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.000160236 + -0.002263921*x + 0.000139687*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.018600433 + 0.004979575*x + 6.62e-06*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.003689806 + 0.000637584*x + -0.000108771*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.020072503 + -0.003865567*x + 3.94e-05*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){0.013722456 + -0.000838848*x + 2.82e-05*x^2}, size=0.5, xlim=c(0,12), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.182109736 + 0.045581451*x + 0.00146424*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.177389815 + 0.038850692*x + 0.001746507*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.138861776 + 0.037070321*x + 0.001705684*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.104224703 + 0.036439517*x + 0.001639454*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.13290001 + 0.046077374*x + 0.00118235*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.182133304 + 0.044615214*x + 0.001391112*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.076156359 + 0.04159499*x + 0.001286267*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.109072887 + 0.04023526*x + 0.00143905*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.03682103 + 0.008410305*x + -0.000205943*x^2}, size=0.5, xlim=c(0,23), colour='#bdbdbd') +
  stat_function(fun=function(x){0.013009709 + -0.028101112*x + 0.001002301*x^2}, size=0.5, xlim=c(0,22), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.017823509 + -0.016588748*x + 0.001951543*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.032498213 + -0.010001844*x + 0.001744399*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046579883 + 0.009738836*x + -0.00123328*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.047868003 + 0.013149289*x + -0.001469638*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.089384041 + -0.005864848*x + 3.97e-05*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.056999468 + -0.01674407*x + 0.000791166*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.066203766 + -0.015095089*x + 0.000525694*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.01504824 + -0.022216259*x + 0.001066483*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025387558 + -0.027398744*x + 0.001303421*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.050425753 + -0.017869209*x + 0.000788277*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.057926773 + -0.01793887*x + 0.000846591*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.069577326 + -0.010052017*x + 0.000399116*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.035896837 + -0.017541434*x + 0.000772584*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.042974859 + -0.00850663*x + 0.000494566*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.062234981 + -0.002058128*x + 2.19e-05*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){0.057485465 + -0.004084372*x + 6.87e-05*x^2}, size=0.5, xlim=c(0,21), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.059752755 + -0.005943841*x + 0.000543005*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.11228537 + -0.00790511*x + 0.000672972*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.011181703 + -0.004278104*x + 0.000463698*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.008951155 + -0.002602959*x + 0.000349383*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.027026976 + -0.005938221*x + 0.000501144*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009679251 + -0.003972973*x + 0.000361442*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.05246031 + -0.014709754*x + 0.000989837*x^2}, size=0.5, xlim=c(0,18), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00038144 + -0.002175523*x + -6.82e-05*x^2}, size=0.5, xlim=c(0,18), colour='#bdbdbd') +
  stat_function(fun=function(x){0.015034417 + 0.000120247*x + 0.000176039*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){0.005414322 + -0.002718822*x + 0.00042538*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025436733 + 0.014442848*x + -5.36e-05*x^2}, size=0.5, xlim=c(0,11), colour='#bdbdbd') +
  stat_function(fun=function(x){0.041540165 + -0.001236289*x + 0.000611532*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017629602 + 0.006626589*x + 0.00018929*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.004196159 + -0.008305614*x + 0.000795488*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.003230513 + -0.0042375*x + 0.000739344*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016881232 + 0.005066636*x + -0.001444389*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.100564279 + 0.029826643*x + -0.002676858*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.105979452 + 0.012976427*x + -0.00194226*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.025670977 + -0.008609304*x + -0.001309907*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.011790404 + 0.009765032*x + -0.002048706*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.028660728 + -0.003193579*x + -0.001676102*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.007301707 + -0.013515089*x + 0.001269964*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.047812815 + -0.022503923*x + 0.001805503*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.015778769 + -0.013199688*x + 0.001548151*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029164492 + -0.003890318*x + 0.001117258*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.079037925 + -0.002541791*x + 0.00063516*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.0059648 + -0.004102002*x + 0.000908388*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.011271639 + -0.00883368*x + 0.001129146*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.061455428 + -0.003396716*x + 0.000724853*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024120603 + -0.00036483*x + 0.000356021*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026408183 + 0.004863091*x + 2.66e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.006595599 + 0.005515608*x + 4.61e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017405209 + -0.023347426*x + 0.003466966*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.074649133 + -0.03406596*x + 0.003530536*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.009132959 + -0.030697748*x + 0.00366252*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.008250086 + -0.010944015*x + 0.005222487*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.015254066 + -0.01776057*x + 0.005440173*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016268422 + -0.025023555*x + 0.005799574*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.001693221 + -0.037724148*x + 0.006245197*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.03734395 + -0.019054265*x + 0.003203527*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.035384281 + -0.025677156*x + 0.003472826*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046574671 + -0.016795194*x + 0.003191933*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.014858307 + -0.018002776*x + 0.003359084*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.015894644 + -0.020252193*x + 0.003361754*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.001527721 + -0.016060148*x + 0.003299181*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046419622 + -0.01282872*x + 0.003028268*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.000992315 + -0.010245772*x + 0.003038629*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009105691 + -0.01787492*x + 0.003426049*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.038857234 + -0.004821049*x + 0.002989843*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.088326659 + -0.008013196*x + 0.002945228*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.078624037 + -0.011221754*x + 0.003223161*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.079755068 + -0.011497543*x + 0.003121896*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.084559788 + -0.016722496*x + 0.003319893*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.116166667 + -0.007018092*x + 0.002916833*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.026439546 + -0.006756865*x + -0.001922767*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.014587352 + -0.017683839*x + -0.001441062*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.032000402 + -0.007517998*x + -0.001633304*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.050088656 + -0.004527321*x + 0.00300267*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.003150274 + 0.006698072*x + 0.002377777*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.00365224 + -0.000442081*x + 0.00208873*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.101463031 + 0.012121865*x + 0.001107273*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.149057904 + 0.020007223*x + 0.000624934*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.146392921 + 0.009754109*x + 0.000798909*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.029980771 + 0.015438137*x + 0.001115005*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.009165033 + 0.023209645*x + 0.000818871*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.032252065 + 0.009780215*x + 0.001327261*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.054032287 + 0.008228787*x + 0.001549411*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.04480628 + 0.004303925*x + 0.001354437*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.072273688 + -0.000245742*x + 0.001729982*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.021251592 + 0.00539538*x + 0.001237236*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.030507362 + 0.000708204*x + -0.000283409*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.048003655 + 0.008276997*x + -0.000716122*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.073751832 + 0.020799331*x + -0.001305697*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.041356607 + 0.01246301*x + -0.000941509*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.101425413 + 0.021188019*x + -0.001387468*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.006098153 + 0.004800486*x + -0.000535199*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.016473553 + 0.00245202*x + 0.00049254*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.01722922 + 0.01510089*x + -0.000109239*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){0.001219399 + 0.006845194*x + -0.000676404*x^2}, size=0.5, xlim=c(0,10), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.021892276 + -0.007545518*x + -6.18e-05*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.002642722 + -0.006618358*x + -0.000189115*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.020483536 + 0.001721552*x + -0.000456843*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.029106448 + -0.007315464*x + 5.18e-05*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.013393065 + 0.005179127*x + -0.001276912*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.001696851 + 0.004359459*x + -0.001194973*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.027453343 + 0.003659475*x + -0.001184372*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.062099073 + 0.004828355*x + -0.001226268*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.061477606 + 0.003056935*x + -0.001179767*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.011933035 + 0.000361858*x + -0.001081095*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){0.024493338 + 0.009008377*x + -0.001446184*x^2}, size=0.5, xlim=c(0,7), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.014574405 + 0.006226715*x + -0.000213818*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.011729989 + 0.007382791*x + -0.000277341*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.017074624 + 0.009238834*x + -7.64e-05*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.01017463 + 0.006196507*x + 3.79e-05*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009177981 + 0.004583857*x + 0.000117192*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029481117 + 0.004154128*x + 0.000179335*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.041444352 + 0.004851058*x + 6.71e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.014649233 + 0.000816602*x + 0.000443105*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.029051003 + 0.000231884*x + 0.00040511*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.005893978 + -0.027479827*x + 0.002832464*x^2}, size=0.5, xlim=c(0,5), colour='#bdbdbd') +
  stat_function(fun=function(x){0.028850903 + -0.024488338*x + 0.002643486*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.000968371 + -0.007331141*x + 0.000474647*x^2}, size=0.5, xlim=c(0,3), colour='#bdbdbd') +
  stat_function(fun=function(x){0.084593553 + -0.001611358*x + -2.08e-05*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.038785035 + -0.016243252*x + 0.000715714*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.030663272 + -0.007359563*x + 0.000254553*x^2}, size=0.5, xlim=c(0,9), colour='#bdbdbd') +
  stat_function(fun=function(x){0.046179229 + -0.007734769*x + 0.000133632*x^2}, size=0.5, xlim=c(0,16), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.026584867 + -0.01217027*x + 0.000738937*x^2}, size=0.5, xlim=c(0,16), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.028051695 + 0.015141046*x + -0.000306*x^2}, size=0.5, xlim=c(0,16), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.014106945 + -0.003892552*x + 0.000527392*x^2}, size=0.5, xlim=c(0,8), colour='#bdbdbd') +
  stat_function(fun=function(x){0.008218949 + -0.000345288*x + 0.00119602*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.042078283 + -0.002065892*x + 0.001051503*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.025795056 + 0.002153583*x + 0.001094266*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){0.047775364 + 0.005394338*x + 0.000831622*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.068824756 + 0.011088776*x + 0.000231654*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.021561419 + 0.014821747*x + -0.000172585*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.052934749 + 0.000796567*x + 0.000649443*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.09633173 + 0.006669682*x + 0.00049931*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.005187837 + 0.022549336*x + -0.000923856*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.072243956 + 0.015875412*x + -0.000340652*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.010021513 + 0.020341071*x + -0.000829221*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.00836217 + -0.002404164*x + -0.000201077*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.032822921 + -0.003077675*x + -0.00016272*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.037871915 + -0.001018463*x + -0.000315469*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.008411798 + -0.003933954*x + -0.000161132*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.02892351 + 0.002044971*x + 0.001221654*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.058116585 + 0.003690254*x + 0.001162055*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.007228781 + 0.001456354*x + -0.001517545*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.020997814 + 0.015313182*x + -0.002257762*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){0.009484363 + 0.011587328*x + -0.002074558*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.044456841 + 0.011674982*x + -0.001912062*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.039234879 + 0.003248388*x + -0.001836572*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.066933796 + 0.007260864*x + -0.001637591*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.056481307 + 0.004141828*x + -0.001821529*x^2}, size=0.5, xlim=c(0,4), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.000199192 + 0.007519739*x + 2.54e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){0.056221637 + 0.009895426*x + -0.000627261*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.020595589 + 0.004910952*x + -3.76e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.104099159 + -0.00253671*x + 0.00052466*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.053806444 + 0.00741595*x + -0.000173121*x^2}, size=0.5, xlim=c(0,6), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.070187613 + 0.003817923*x + 0.000153145*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.031836741 + 0.003627604*x + 4.82e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.060027139 + 0.004714096*x + 4.29e-05*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.126135123 + -0.000588978*x + 0.00052623*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.079409497 + 0.010291094*x + -0.000190674*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.227560838 + 0.00392271*x + 0.000586313*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.131284251 + 0.005102001*x + 0.000164646*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
  stat_function(fun=function(x){-0.064842918 + 0.011684862*x + -0.000238832*x^2}, size=0.5, xlim=c(0,2), colour='#bdbdbd') +
#estimated as mean across treatment lines
  #overall line (because plot mani not significant)
  stat_function(fun=function(x){(-0.02694475 + -0.002949595*x + 0.00675485*x^2)*0.09064568 - 0.00235573}, size=3, xlim=c(0,23), colour='black')

  

dispersionInsetPlot <- function() {
  print(dispersionPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(dispersionInset, vp=viewport(width = 0.13, height = 0.25, x = 0.83, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

dispersionInsetPlot() #export at 1200x1000







###richness
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
richnessInset <- ggplot(data=rawData, aes(x=S_PC)) +
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
stat_function(fun=function(x){-0.021677338 + -0.007322007*x + 0.000567061*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.018172345 + 0.006508211*x + -0.000431271*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.096780009 + -0.024208862*x + 0.002067667*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.017835229 + -0.033655242*x + 0.00051909*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.066154966 + -0.027990936*x + 0.003885528*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){-0.045445831 + -0.021292455*x + 0.003452466*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.010734327 + -0.016600072*x + 0.002688641*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.000309782 + -0.022942987*x + 0.003128759*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.07652313 + -0.032779878*x + -0.000376624*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.081341003 + 0.013908206*x + -0.002606891*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0.023000503 + -0.051883616*x + 0.001326693*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.062401703 + -0.019986831*x + -0.001075759*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.116508313 + -0.068095166*x + 0.002493991*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.159103355 + -0.089096458*x + 0.003843*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.20294003 + -0.064663352*x + -0.002114603*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.169927674 + -0.029430969*x + -0.004460169*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.190826393 + -0.041269132*x + -0.003500128*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.20369155 + -0.088901145*x + 0.011081314*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.0341075 + -0.046223084*x + 0.00793682*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.038339548 + -0.052007802*x + 0.00848739*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.012831191 + -0.004681017*x + 0.005802431*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.008659413 + -0.039544248*x + 0.007704011*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.149098024 + -0.077034509*x + 0.01058593*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.161159402 + -0.072990112*x + 0.010050752*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.192147157 + -0.074718998*x + 0.010303687*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.041279427 + 0.023273235*x + -0.00085151*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.057459101 + 0.002152223*x + 0.000833243*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011299454 + 0.015148455*x + -0.000194803*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.065861648 + 0.021935787*x + -0.000935321*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.001618991 + 0.012187279*x + -1.39e-05*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.02295008 + -0.016294204*x + 0.000427408*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.033512305 + -0.005615157*x + -0.000434468*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.042908483 + -0.014813072*x + 2.6e-05*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.045069962 + -0.031694564*x + 0.001287275*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){-0.010170109 + -0.054695185*x + 0.003253768*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.048015667 + -0.027549195*x + 0.001426462*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.136630589 + -0.010103674*x + 0.000167738*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.020796263 + 0.027382012*x + -0.005741985*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.034824372 + 0.04016861*x + -0.005971992*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.078674491 + 0.007515822*x + -0.004118268*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.017963189 + 0.028087563*x + -0.005538724*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.008567582 + 0.031221262*x + -0.00556602*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.046603421 + 0.053060982*x + -0.006878391*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.07964425 + -0.005831557*x + 0.001281166*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.009718768 + -0.012599049*x + 0.002003*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.073554858 + -0.033597494*x + 0.003569036*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.010388979 + -0.032316639*x + 0.002962308*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.032178841 + -0.021840143*x + 0.002761689*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.034723742 + -0.016412604*x + 0.002887236*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.019965839 + 0.045756757*x + -0.000755811*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.07200993 + -0.018359524*x + 0.00318964*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.069409088 + 0.003210692*x + 0.001844952*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.026959822 + -0.016285753*x + 0.002872976*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.09698491 + -0.028272585*x + 0.001439422*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.122271653 + -0.013343262*x + 0.000497483*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.051309 + -0.036925819*x + 0.002139475*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.057974633 + -0.067575551*x + 0.004237366*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.062389125 + -0.027518205*x + 0.001604813*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.0205894 + -0.018982765*x + 0.001918509*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.024182444 + -0.022026125*x + 0.002383859*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.072108604 + -0.042744042*x + 0.00375886*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.083292938 + -0.054479746*x + 0.004263495*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.017285082 + -0.028799861*x + 0.002722009*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.10363195 + -0.010676588*x + 0.00121743*x^2}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.115565253 + -0.019964087*x + 0.001978936*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){-0.053708554 + -0.008573657*x + 0.00111012*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){-0.003441759 + -0.050991785*x + 0.002516999*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){-0.328470365 + -0.083460855*x + 0.004448551*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.12349762 + 0.011226895*x + -0.000850106*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){-0.127092836 + -0.076552058*x + 0.002403951*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.031874443 + -0.010249449*x + 0.000235037*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){-0.282024181 + -0.073088569*x + 0.003265135*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.000445975 + -0.013622361*x + 0.000730122*x^2}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){-0.286311118 + -0.063842534*x + 0.002679653*x^2}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){-0.005347891 + -0.050118137*x + 0.003486886*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){-0.331620759 + -0.123083434*x + 0.00705288*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.013641266 + -0.021222563*x + 0.00193886*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){-0.209323951 + -0.121476219*x + 0.006384082*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){-0.04194815 + -0.022822678*x + 0.002679459*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){-0.305382265 + -0.097166384*x + 0.005974897*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.09024787 + 0.006223052*x + -0.003226072*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.06409514 + 0.006514637*x + -0.003291058*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.028767503 + -0.002542442*x + -0.00091492*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.034525119 + -0.01073126*x + -6.67e-05*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){-0.087888385 + -0.044584879*x + 0.001122572*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.056062155 + -0.044974133*x + 0.000735419*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.067610162 + -0.049142144*x + 0.005701355*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.196543078 + -0.005866948*x + 0.002014397*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.088434821 + -0.022976264*x + 0.003192418*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.038273002 + -0.034161592*x + 0.001081905*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){-0.008161862 + -0.089926767*x + 0.003903584*x^2}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){-0.21376675 + -0.068819585*x + 0.003280229*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.188736841 + -0.069380366*x + 0.003020113*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.062438067 + 0.015676418*x + -0.002547737*x^2}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.066811277 + -0.092173896*x + 0.0042668*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.122903448 + -0.033134941*x + 0.003640883*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.008145967 + -0.03803995*x + 0.004555676*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.026624549 + -0.014879008*x + 0.003416342*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.126663337 + -0.020205609*x + 0.002347026*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.009614823 + 0.001410203*x + 0.001858533*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.042036322 + 0.00908473*x + 0.001051133*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011881277 + -0.004809594*x + 0.002150447*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.021080764 + -0.011509561*x + 0.002757435*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.089913963 + 0.017097011*x + 0.000382858*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.029596548 + 0.02623026*x + -0.002507142*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){-0.030153212 + 0.015307519*x + -0.001913016*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.01397403 + 0.001418953*x + -0.001230633*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.070134863 + 0.020158885*x + -0.00261828*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.037242463 + 0.005291953*x + -0.001823333*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){-0.029509297 + 0.018707886*x + -0.002298656*x^2}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){-0.061745443 + 0.024146655*x + -0.002120516*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.040619839 + 0.011396547*x + -0.001407808*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.033891953 + 0.009796628*x + -0.002176754*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.027543709 + 0.022690327*x + -0.002711683*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.009819284 + -0.005170591*x + -0.000889404*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.022951789 + -0.0001754*x + -0.000795599*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){-0.029697863 + 0.02868894*x + -0.002484557*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){-0.026841697 + 0.010413933*x + -0.0015047*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.055875871 + 0.041780628*x + -0.002163004*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.094887811 + 0.033769593*x + -0.001839344*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.043857604 + 0.032008803*x + -0.001117424*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.023367579 + 0.03110691*x + -0.00100016*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.074068224 + 0.032359749*x + -0.001568856*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.080227902 + 0.030756307*x + -0.001444645*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.019747421 + 0.020770576*x + -0.000591734*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.061064467 + 0.005884993*x + 0.00044606*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.045932633 + -0.010018196*x + -0.000460719*x^2}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){-0.130744892 + -0.039476895*x + 0.00168356*x^2}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){0.003333246 + -0.031452253*x + 0.004394314*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.007003369 + -0.060377787*x + 0.005546626*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.00899515 + 0.046889298*x + -0.003497132*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.121173869 + 0.00534428*x + -0.000346315*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.078959977 + -0.061906785*x + 0.003076305*x^2}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){-0.100379645 + -0.044753662*x + 0.002846664*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0.066812466 + -0.038452188*x + 0.002755388*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0.137240993 + -0.08940921*x + 0.005854724*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0.145142822 + -0.061484027*x + 0.004016289*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.121665811 + -0.04013053*x + 0.002592162*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.035270344 + -0.009302673*x + 0.000349779*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){-0.094697364 + -0.025218704*x + 0.00122207*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.065711406 + -0.003205008*x + 0.000578554*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.066747354 + -0.026503676*x + 0.000810205*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.057863175 + -0.013928023*x + 0.000243239*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){-0.028107686 + -0.009385395*x + 0.000661807*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){-0.010764967 + 0.031110912*x + -0.002517045*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.176030405 + -0.020277342*x + 0.001342659*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.117342402 + 0.016998577*x + -0.001822358*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011234731 + -0.005974919*x + 0.000195927*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.002052384 + 0.012118496*x + -0.001427571*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.03918175 + -0.009565465*x + 0.00013892*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.094577329 + -0.000771254*x + -0.000218829*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){-0.073809794 + -0.02085897*x + 0.000460574*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){-0.051285323 + 0.00417876*x + -0.000599421*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0.04424445 + -0.003670947*x + -0.000371596*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0.082894764 + -0.012229818*x + -0.000127118*x^2}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){0.090323114 + 0.033779885*x + -0.003453232*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.094726918 + 0.01882863*x + -0.002339971*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.121733965 + 0.002360852*x + -0.000317063*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.085031772 + 0.013890127*x + -0.001289006*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.013808106 + -0.001321949*x + -0.001996972*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.076164861 + -0.032952435*x + 0.000171806*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.057883403 + -0.058565996*x + 0.001067463*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.043315462 + 0.025278853*x + -0.004092482*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.048546596 + -0.028524731*x + -0.000967532*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.017160439 + -0.038680434*x + -0.000912038*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.085306868 + -0.018014719*x + -0.000531566*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.136783014 + -0.061181566*x + 0.002438748*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.023703545 + -0.007775161*x + -0.000580209*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.024852203 + -0.031836132*x + 0.001346202*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.144074095 + -0.007490402*x + 0.001995401*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.160032813 + -0.007054081*x + 0.001869775*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.124682191 + -0.014709127*x + 0.002489372*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.160997943 + -0.004874866*x + 0.001756558*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.089702183 + -0.054107874*x + 0.002548228*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.145877767 + -0.062459448*x + 0.002253471*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.069240288 + -0.031214973*x + 0.00129165*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.152396504 + -0.02019422*x + 0.007380956*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0.070803152 + -0.040973075*x + 0.007602615*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.130729798 + -0.041330082*x + 0.008157942*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.234994002 + -0.065066671*x + 0.008411243*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0.181823928 + -0.06072736*x + 0.008017049*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.173282301 + -0.039074491*x + 0.007135362*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.079904257 + -0.037526167*x + 0.005962066*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.120389987 + -0.062374256*x + 0.009246229*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.056066839 + -0.043190015*x + 0.007183138*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.097411574 + -0.026831408*x + 0.007177524*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.036501456 + -0.044130788*x + 0.007197329*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.100884303 + -0.021056456*x + 0.006479669*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.052403913 + -0.044464238*x + 0.007257364*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.019223003 + -0.056466152*x + 0.008681354*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.035463948 + -0.064391309*x + 0.008942751*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.039765294 + -0.05349369*x + 0.0080208*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.103164822 + -0.061970479*x + 0.009297573*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){-0.065478357 + -0.046604104*x + 0.008202677*x^2}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){-0.013926568 + -0.028987741*x + 0.006762392*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.01897864 + -0.045453839*x + 0.007324294*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.027172866 + -0.036390653*x + 0.006757898*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.164054568 + -0.032235677*x + 0.005631738*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.293921699 + 0.098790523*x + -0.000214564*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.119439043 + 0.023438016*x + 0.004802367*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.292126374 + 0.103761056*x + -0.000619236*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.155884807 + -0.021774162*x + 0.001814247*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.125228901 + -0.018263732*x + 0.001449515*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.047006848 + -0.0439309*x + -0.00270218*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.018066684 + -0.071291643*x + -0.002095642*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.023848428 + -0.080299254*x + -0.001521272*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.091251599 + -0.101180477*x + -0.002004174*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.082439353 + -0.0378053*x + -0.003712491*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.074040545 + -0.063859687*x + -0.00195594*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.077165054 + -0.054361964*x + -0.002676154*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.160746369 + 0.003755281*x + -0.006900301*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.102309699 + -0.042242724*x + -0.0032494*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.098693904 + -0.038721487*x + -0.00353592*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.020050911 + -0.038785752*x + -0.003722199*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.076284135 + 0.003683274*x + -0.001162814*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.129185361 + -0.041465245*x + 0.002045595*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){-0.080326513 + -0.027044217*x + 0.001006764*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){-0.109489167 + -0.048692627*x + 0.00247231*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){-0.119577632 + -0.048443225*x + 0.002385414*x^2}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){-0.018964111 + 0.009103346*x + -0.001107093*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.028350461 + 0.012508814*x + -0.001374468*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.156562234 + 0.013930105*x + -0.003569162*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.044519015 + 0.024257348*x + -0.002740088*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.017182286 + -0.012629524*x + -0.000626615*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.029049151 + -0.03147215*x + 0.000349889*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.11045396 + 0.009246401*x + -0.002493054*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.068894848 + -0.028863212*x + 0.000529979*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){-0.054555673 + -0.014974229*x + -0.001006945*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.063627327 + 0.023438245*x + -0.003590923*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.020510262 + -0.001042093*x + -0.001637787*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.001849181 + -0.001591569*x + -0.002106212*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.018174969 + -0.001021019*x + -0.0020155*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.031505773 + 0.043652454*x + -0.004862013*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){-0.025548606 + 0.00415095*x + -0.00205998*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){0.088315667 + -0.00259055*x + -0.00081508*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.097580112 + 0.00586298*x + -0.001316569*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.058520258 + -0.014323909*x + 0.000423657*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.075802942 + 0.002319536*x + -0.000880516*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.093180424 + 0.006216385*x + -0.001154659*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.055388257 + -0.026629348*x + 0.001277353*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.05075954 + -0.003140308*x + -0.000594801*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.083531464 + -0.03065705*x + 0.001624594*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.035191333 + -0.01764837*x + 0.000672772*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.139736951 + -0.057179022*x + 0.006835795*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.026462969 + -0.034745815*x + 0.004697403*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.120191015 + -0.055673808*x + 0.003310635*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.168289018 + 0.029145776*x + -0.002974509*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){-0.033962845 + -0.037702383*x + 0.001956475*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.116965155 + 0.011113595*x + -0.00132176*x^2}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){0.050216597 + 0.002032879*x + -0.000236908*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.057770567 + 0.012972625*x + -0.000758053*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0.000909437 + 0.012012961*x + 4.23e-08*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0.022955448 + -0.003498435*x + 0.000555244*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.051123172 + -0.045183282*x + 0.001148106*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.073442324 + -0.016891109*x + -9.38e-05*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.017493317 + -0.041115786*x + 0.001109189*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.078386837 + -0.018368181*x + -0.00099763*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.06496016 + -0.031212343*x + -0.000207085*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.032259923 + -0.007202144*x + -0.001762261*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.059198022 + -0.02497948*x + -0.000840215*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.101819479 + -0.024414124*x + -0.001417238*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.096039907 + -0.005355727*x + -0.00176928*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.116357375 + -0.010027938*x + -0.001815857*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.035349018 + -0.027886418*x + 0.000142491*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.043826608 + -0.010152376*x + 0.000248084*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.060611896 + -0.028408206*x + 0.001553407*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.038173395 + -0.045897981*x + 0.00257732*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.024395488 + 0.005616415*x + -0.000917709*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.016728074 + -0.01925411*x + 0.001896801*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.079791421 + 0.01989878*x + -0.000989734*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.068256748 + 0.007025001*x + -0.00190558*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.115979442 + -0.026617799*x + -0.000755203*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.073630943 + -0.014727331*x + -0.00119967*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.083969266 + -0.01787497*x + -0.001129681*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.03844796 + -0.039719321*x + -0.0018645*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.197549374 + -0.006398479*x + -0.002283233*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){0.009026597 + -0.034835696*x + -0.001962263*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.003165492 + 0.046179058*x + -0.005080928*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.110014848 + 0.018963165*x + -0.001168549*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.041992639 + 0.034363994*x + -0.001200398*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.028205006 + 0.01316963*x + 6.63e-05*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.019189047 + -0.015506594*x + 0.001259998*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.042602249 + 0.012431329*x + 0.000110983*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.007389078 + 0.011436799*x + 0.000277376*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.007232073 + 0.016654218*x + -0.000368221*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.120287757 + -0.021413096*x + 0.002805772*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.005953613 + 0.01214737*x + 0.000185074*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){-0.166418787 + -0.005333977*x + 0.00184107*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.071263876 + -0.001682708*x + 0.00131175*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){-0.300856218 + -0.029632228*x + 0.003821639*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(0.319741 + -0.0520446*x + 0.00197122*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(0.319741 + -0.0520446*x + 0.00197122*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(0.319741 + -0.0520446*x + 0.00197122*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(0.319741 + -0.0520446*x + 0.00197122*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((0.319741-0.694714) + (-0.0520446-0.3041105)*x + (0.00197122+0.009461835)*x^2)*0.2287037 - 0.07758351}, size=3, xlim=c(0,22), colour='#EC1804')
                                  


richnessInsetPlot <- function() {
  print(richnessPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(richnessInset, vp=viewport(width = 0.13, height = 0.25, x = 0.30, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

richnessInsetPlot() #export at 1200x1000






###evenness
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
evennessInset <- ggplot(data=rawData, aes(x=SimpEven_change)) +
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
  stat_function(fun=function(x){0.007564893 + 0.02967549*x + -0.001713149*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.002372168 + 0.025216304*x + -0.001405815*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.086278093 + 0.017750593*x + -0.001452646*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.082879172 + 0.008724295*x + -0.001448261*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.05019695 + -0.002997139*x + 5.36e-05*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.017609879 + -0.011494803*x + 0.000685837*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.046372123 + -0.003018496*x + 1.5e-05*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.028070941 + -0.009600944*x + 0.000465306*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.005934857 + -0.012809289*x + 0.000268991*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.029900955 + 0.004388769*x + -0.000635979*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){-0.043284622 + -0.031154267*x + 0.001396124*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.022936805 + 0.002777701*x + -0.000675217*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.014987202 + -0.003046121*x + -0.000423464*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.032099107 + 0.00408971*x + -0.000897488*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.009041208 + 0.015086716*x + -0.002430533*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.007665843 + 0.013395711*x + -0.002329999*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.005715601 + 0.010894626*x + -0.002120862*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008439051 + -0.005768154*x + -0.000366219*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.015352128 + -0.003753727*x + -0.000622794*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.01190844 + -0.000331102*x + -0.000714776*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.008207993 + -0.000665532*x + -0.000637974*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.033165423 + 0.005597321*x + -0.001325169*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.001362937 + -0.005042951*x + -0.000558413*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.016833005 + -0.005840469*x + -0.000306527*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.0112123 + -0.004569071*x + -0.00041819*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.046155757 + -0.000427218*x + 6.62e-05*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.034409424 + -0.002987117*x + 0.000273401*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.026016829 + -0.004698203*x + 0.000443266*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.047001726 + -0.000682805*x + 7.08e-05*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.031568534 + -0.005268165*x + 0.000446189*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.036519302 + 0.001128754*x + -0.000472404*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.025855904 + -0.001810984*x + -0.000330515*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.00549804 + -0.00801669*x + 0.000130336*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.009044503 + -0.021786924*x + 0.001164979*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){-0.064355067 + -0.031615234*x + 0.001953928*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.028336341 + -0.011992486*x + 0.00066197*x^2}, size=0.5, xlim=c(0,15), colour='#4A06AC44') +
  stat_function(fun=function(x){0.000301999 + 0.000100299*x + -0.000136238*x^2}, size=0.5, xlim=c(0,15), colour='#1400E544') +
  stat_function(fun=function(x){0.004706577 + -0.012840988*x + 0.002229226*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.007943385 + -0.015043846*x + 0.002526124*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.03841933 + -0.002326306*x + 0.001765482*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.00931768 + -0.010019781*x + 0.00209756*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.022541209 + -0.025946022*x + 0.003086287*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.010169585 + -0.013275375*x + 0.002288334*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.038968628 + 0.011550912*x + -0.001770808*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.049013614 + 0.00664273*x + -0.001458624*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.072974282 + 0.004140357*x + -0.001155096*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.049905348 + 0.011499509*x + -0.001794348*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.067936379 + 0.00088109*x + -0.001015079*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.03287821 + 0.004564091*x + -0.001330982*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.012395724 + 0.003951218*x + -0.001163965*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.016877007 + 0.006331455*x + -0.00132096*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.008786451 + 0.004243245*x + -0.001265012*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.025142703 + 0.004426995*x + -0.001271145*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.036853625 + 0.007700539*x + -0.001363592*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.024900657 + 0.009235408*x + -0.001376054*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.023481226 + 0.005560388*x + -0.001175952*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.011669491 + 0.009993035*x + -0.001411509*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.025395237 + 0.010766792*x + -0.001445898*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.002029003 + 0.013417637*x + -0.001848129*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.017877861 + 0.010346523*x + -0.001596785*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.009935657 + 0.017021237*x + -0.002126633*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.023321728 + 0.018973598*x + -0.002400245*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.014169956 + 0.008889791*x + -0.001535764*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.005299153 + 0.022042385*x + -0.00170563*x^2}, size=0.5, xlim=c(0,13), colour='#4A06AC44') +
  stat_function(fun=function(x){0.013719672 + 0.021564663*x + -0.001696978*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){-0.018244659 + 0.012914834*x + -0.000809716*x^2}, size=0.5, xlim=c(0,13), colour='#1400E544') +
  stat_function(fun=function(x){0.017112098 + 0.011239187*x + -0.00043966*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.070675543 + 0.05034904*x + -0.003066594*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){-0.014035141 + -0.001356031*x + 0.000148978*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.01848311 + 0.035473984*x + -0.001257482*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.004948489 + -0.000834508*x + 9.88e-05*x^2}, size=0.5, xlim=c(0,22), colour='#B6123C44') +
  stat_function(fun=function(x){0.004376702 + 0.032369929*x + -0.001469313*x^2}, size=0.5, xlim=c(0,22), colour='#EC180444') +
  stat_function(fun=function(x){0.047965237 + 0.005940042*x + -0.000463148*x^2}, size=0.5, xlim=c(0,20), colour='#B6123C44') +
  stat_function(fun=function(x){0.006473859 + 0.017518661*x + -0.000870656*x^2}, size=0.5, xlim=c(0,20), colour='#EC180444') +
  stat_function(fun=function(x){0.010017425 + 0.007141886*x + -0.000635799*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.164597191 + 0.093763562*x + -0.006221038*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){-0.021927586 + 0.009998641*x + -0.00068973*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.011424538 + 0.065170315*x + -0.003456405*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.004199224 + 0.008147414*x + -0.000513364*x^2}, size=0.5, xlim=c(0,9), colour='#B6123C44') +
  stat_function(fun=function(x){0.023189866 + 0.048082421*x + -0.00239693*x^2}, size=0.5, xlim=c(0,9), colour='#EC180444') +
  stat_function(fun=function(x){0.019148564 + 0.001252523*x + 0.000290569*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.024165628 + 0.002740267*x + 7.26e-05*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.030067062 + 0.004007073*x + -0.000392035*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.025454444 + 0.003811976*x + -0.000490427*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){-0.007720244 + 0.016154945*x + -0.00180316*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.004244168 + 0.020472232*x + -0.002132239*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.013916616 + -0.001480002*x + -0.000444885*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.009266558 + -0.006583555*x + -0.000213039*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.0200077 + -0.00300999*x + -0.000424403*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){0.022016174 + 0.002392912*x + -0.001321803*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.016324089 + 0.005936293*x + -0.001321787*x^2}, size=0.5, xlim=c(0,7), colour='#B6123C44') +
  stat_function(fun=function(x){0.036384125 + 0.009460943*x + -0.000901554*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.003358588 + 0.004110488*x + -0.000470901*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.003902289 + 0.006364696*x + -0.000866544*x^2}, size=0.5, xlim=c(0,7), colour='#4A06AC44') +
  stat_function(fun=function(x){0.046735871 + 0.02387679*x + -0.002059546*x^2}, size=0.5, xlim=c(0,7), colour='#800C7444') +
  stat_function(fun=function(x){0.018956434 + -0.006375008*x + 6.31e-05*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.035090221 + -0.013622885*x + 0.000410285*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.01344243 + -0.012032305*x + 0.000370405*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.029867192 + -0.009049439*x + 0.000464487*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.035814619 + 0.006863155*x + -0.000504074*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.03645022 + 0.010063334*x + -0.000764667*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.062646686 + 0.000998915*x + 5.33e-06*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.060674105 + 0.001211417*x + -2.19e-05*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.05696377 + 0.004155959*x + -0.00030076*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.004536649 + -0.00405274*x + 0.000118285*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.013784029 + 0.001360158*x + -0.000302478*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.010753649 + -0.004928422*x + 7.48e-06*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.023084165 + 0.000398282*x + -0.00044783*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.008339343 + -0.002014991*x + -1.52e-05*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.020391009 + -0.002808713*x + -0.000227889*x^2}, size=0.5, xlim=c(0,12), colour='#B6123C44') +
  stat_function(fun=function(x){0.00543614 + -0.000914439*x + -0.000109832*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.013084733 + -0.004915592*x + -7.7e-05*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011011069 + -0.00162443*x + -0.000189416*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.001790136 + -0.005734825*x + 0.000205255*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.025987648 + -0.004160887*x + -0.00035772*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.014659261 + -0.002287434*x + 7.65e-05*x^2}, size=0.5, xlim=c(0,12), colour='#800C7444') +
  stat_function(fun=function(x){0.008079848 + -0.001738023*x + 0.000173413*x^2}, size=0.5, xlim=c(0,12), colour='#1400E544') +
  stat_function(fun=function(x){0.011882609 + -0.004307673*x + 0.000153842*x^2}, size=0.5, xlim=c(0,12), colour='#4A06AC44') +
  stat_function(fun=function(x){0.032393455 + 0.003884886*x + -0.000369428*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.014790027 + 0.001056434*x + -0.000111511*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.056744347 + 0.009996681*x + -0.000989181*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.02755416 + 0.002952692*x + -0.00044792*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.013057889 + -0.001951839*x + 0.00015802*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.047267167 + 0.007993558*x + -0.000602847*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.020108564 + 0.006157354*x + -0.000465138*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.025643939 + 0.006743554*x + -0.000540645*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.023272839 + -0.010403562*x + 0.000516358*x^2}, size=0.5, xlim=c(0,23), colour='#1400E544') +
  stat_function(fun=function(x){0.044718093 + -0.012132111*x + 0.00046048*x^2}, size=0.5, xlim=c(0,22), colour='#4A06AC44') +
  stat_function(fun=function(x){0.026335617 + 0.011040611*x + -0.001018194*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.0415344 + 0.007628348*x + -0.001069481*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.048183884 + -0.008747959*x + 0.000860225*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.046039255 + -0.010499489*x + 0.001048692*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.041167188 + 0.003426117*x + -0.000546895*x^2}, size=0.5, xlim=c(0,21), colour='#B6123C44') +
  stat_function(fun=function(x){0.049862421 + -0.001895043*x + -0.000366548*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.048795743 + -0.006995436*x + -0.000141798*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.037763975 + -0.008299449*x + 6.37e-05*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.011958808 + -0.011563426*x + 0.000280893*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.014533098 + -0.003755015*x + -0.000166219*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.048663255 + -0.00252057*x + -0.000401707*x^2}, size=0.5, xlim=c(0,21), colour='#800C7444') +
  stat_function(fun=function(x){0.017623123 + -0.006885438*x + 0.000363337*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.029332415 + -0.008048059*x + 2.05e-05*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.008376162 + -0.013655774*x + 0.000619645*x^2}, size=0.5, xlim=c(0,21), colour='#4A06AC44') +
  stat_function(fun=function(x){0.008310229 + -0.014781456*x + 0.001129917*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0.022528404 + -0.006834258*x + 5.52e-05*x^2}, size=0.5, xlim=c(0,21), colour='#1400E544') +
  stat_function(fun=function(x){0.021650146 + 0.018631622*x + -0.001071606*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.016244411 + 0.003854737*x + 6.78e-05*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.043102639 + 0.018130836*x + -0.001278183*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.01082511 + 0.007957235*x + -0.000421136*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.039262195 + 0.019056649*x + -0.001132343*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.015051523 + 0.007705483*x + -0.00027624*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.012797035 + 0.014033121*x + -0.000770548*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0.00580361 + 0.006828601*x + -5.39e-05*x^2}, size=0.5, xlim=c(0,18), colour='#1400E544') +
  stat_function(fun=function(x){0.022439044 + 0.007120027*x + -0.000234839*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.035028838 + 0.012033701*x + -0.000747588*x^2}, size=0.5, xlim=c(0,11), colour='#1400E544') +
  stat_function(fun=function(x){0.045413215 + -0.001365593*x + 1.69e-05*x^2}, size=0.5, xlim=c(0,11), colour='#4A06AC44') +
  stat_function(fun=function(x){0.016914483 + -0.004839311*x + 0.000551587*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.001387552 + 0.001001314*x + 0.000240037*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.004820765 + 0.003142483*x + 0.000216352*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.000804342 + -0.000100357*x + 0.000246805*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.009500314 + -0.0064492*x + 0.000146045*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.024372981 + 0.003389794*x + -0.000495484*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.056572351 + 0.018205406*x + -0.001539311*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.015153584 + 0.006342827*x + 0.000472382*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.028610212 + 0.004362246*x + 0.000654714*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){-0.01262416 + 0.017310931*x + -7.6e-05*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.004120103 + 0.001460015*x + -0.000507245*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.020034562 + 0.00221273*x + -0.000711546*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.054011176 + -0.003502936*x + -0.00058058*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.057744264 + -0.007716689*x + -0.000474142*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){0.003862936 + 0.011712772*x + -0.001008825*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.034889892 + 0.004005366*x + -0.00034148*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.04448384 + -0.000245158*x + -7.19e-05*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.030402005 + 0.004254699*x + -0.000387219*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.004786841 + -0.005832841*x + -0.000960065*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){-0.00934816 + -0.00158999*x + -0.00097922*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.015276039 + -0.003833287*x + -0.000820894*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.020789182 + -0.000512896*x + -0.000769652*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0.010672936 + 0.003635662*x + -0.000927174*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.020965833 + -0.002968056*x + -0.000503304*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.030769625 + 0.018556173*x + -0.001466841*x^2}, size=0.5, xlim=c(0,8), colour='#B6123C44') +
  stat_function(fun=function(x){-0.000387124 + 0.012599904*x + -0.000860151*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.004024648 + 0.009000653*x + -0.000788478*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){-0.016012945 + 0.003470988*x + -0.000358971*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.021461038 + 0.009313887*x + -0.001296991*x^2}, size=0.5, xlim=c(0,8), colour='#800C7444') +
  stat_function(fun=function(x){0.025388618 + 0.006809875*x + -0.001293919*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.003283547 + 0.006616558*x + -0.001125362*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.011194298 + 0.000196784*x + -0.000364838*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){0.012637615 + 0.000334785*x + -0.000380226*x^2}, size=0.5, xlim=c(0,8), colour='#4A06AC44') +
  stat_function(fun=function(x){0.01846888 + -0.001095725*x + -0.000392018*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.002857538 + -0.013392556*x + 2.31e-05*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.067314062 + -0.001839824*x + -0.000979106*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.007949274 + -0.014361756*x + -1.92e-05*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.011229727 + -0.01070484*x + 9.04e-06*x^2}, size=0.5, xlim=c(0,3), colour='#800C7444') +
  stat_function(fun=function(x){0.010771981 + -0.004960919*x + -0.000454291*x^2}, size=0.5, xlim=c(0,3), colour='#B6123C44') +
  stat_function(fun=function(x){0.004049309 + -0.003530494*x + -0.000595634*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.04111651 + -0.000465534*x + -0.000373831*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.042102929 + 0.004517022*x + -0.000667754*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.022819314 + 0.001414032*x + -0.000477037*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.12145939 + 0.020537545*x + -0.000478485*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.118716031 + 0.014361601*x + -0.000115718*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.115461233 + 0.013045564*x + -0.000142664*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.000378224 + -0.021438277*x + 0.001422074*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.000968576 + -0.020463286*x + 0.001381086*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.015558649 + 0.030534594*x + -0.001167921*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.041383502 + 0.054101107*x + -0.002664611*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.058818337 + 0.059377975*x + -0.003073234*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.060921803 + 0.089161183*x + -0.004407644*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.009822883 + 0.03309308*x + -0.001419524*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.006291365 + 0.037871127*x + -0.001585837*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.012493586 + 0.036455543*x + -0.001489279*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008640162 + 0.03457108*x + -0.001325541*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.005193458 + 0.037705698*x + -0.001405996*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.007391041 + 0.04268077*x + -0.001782785*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.014449457 + 0.044984623*x + -0.001941207*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.019038085 + 0.016010719*x + -0.001256386*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.017091873 + 0.023346373*x + -0.001558192*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){0.00882989 + 0.024587064*x + -0.001550011*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0.025174486 + 0.025985582*x + -0.00170006*x^2}, size=0.5, xlim=c(0,5), colour='#B6123C44') +
  stat_function(fun=function(x){0.000200229 + 0.021055758*x + -0.001329353*x^2}, size=0.5, xlim=c(0,5), colour='#EC180444') +
  stat_function(fun=function(x){0.000369795 + 0.016266231*x + -0.000981826*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.011051933 + 0.007053453*x + 0.001635781*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.017553197 + 0.048633988*x + -0.002858259*x^2}, size=0.5, xlim=c(0,10), colour='#4A06AC44') +
  stat_function(fun=function(x){0.04224757 + -0.010261776*x + 0.00170141*x^2}, size=0.5, xlim=c(0,10), colour='#1400E544') +
  stat_function(fun=function(x){0.023702458 + -0.027548495*x + 0.002674845*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.009414409 + -0.03411544*x + 0.003109371*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.033825173 + -0.021543926*x + 0.002182266*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.009431107 + -0.033002789*x + 0.003116068*x^2}, size=0.5, xlim=c(0,5), colour='#800C7444') +
  stat_function(fun=function(x){-0.022996074 + -0.001560157*x + -0.000267215*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008113209 + 0.004316868*x + -0.000759185*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.007312179 + 0.003801468*x + -0.000631524*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.051524835 + -0.005761974*x + 0.000345076*x^2}, size=0.5, xlim=c(0,5), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.049843448 + -0.007171502*x + 0.000397599*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){-0.043106886 + 0.003764923*x + -0.000330558*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){-0.030977878 + -0.001157785*x + -0.000251536*x^2}, size=0.5, xlim=c(0,7), colour='#1400E544') +
  stat_function(fun=function(x){-0.012478175 + 0.013712658*x + -0.001519867*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.008902148 + 0.018334326*x + -0.001808564*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.021807891 + 0.005422517*x + -0.000167572*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){-0.005065056 + 0.011853896*x + -0.000426285*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.000458999 + 0.013938435*x + -0.000595008*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.036546302 + 0.007020553*x + -7.22e-05*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.022600678 + 0.010778686*x + -0.00035051*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-7.25e-05 + 0.013712348*x + -0.000620643*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){-0.012468505 + 0.009323439*x + -0.00044728*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.008068891 + 0.009867358*x + -0.000557353*x^2}, size=0.5, xlim=c(0,5), colour='#1400E544') +
  stat_function(fun=function(x){0.01470014 + 0.009113981*x + -0.000649231*x^2}, size=0.5, xlim=c(0,3), colour='#4A06AC44') +
  stat_function(fun=function(x){0.012395289 + 0.005900145*x + -0.000410686*x^2}, size=0.5, xlim=c(0,3), colour='#1400E544') +
  stat_function(fun=function(x){0.023433738 + 0.010823055*x + -0.000884698*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){0.01821192 + 0.022532363*x + -0.001136707*x^2}, size=0.5, xlim=c(0,9), colour='#1400E544') +
  stat_function(fun=function(x){-0.032407889 + -0.002513289*x + 0.000345662*x^2}, size=0.5, xlim=c(0,9), colour='#4A06AC44') +
  stat_function(fun=function(x){0.004031778 + 0.009019021*x + -0.000158562*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){0.07658946 + 0.021887919*x + -0.00146121*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0.010367987 + 0.011320014*x + -0.000805368*x^2}, size=0.5, xlim=c(0,16), colour='#1400E544') +
  stat_function(fun=function(x){-0.00423093 + 0.009392951*x + -0.00081311*x^2}, size=0.5, xlim=c(0,8), colour='#1400E544') +
  stat_function(fun=function(x){-0.015245477 + 0.017808097*x + -0.001890435*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.023020482 + 0.012028427*x + -0.001438751*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){-0.006082813 + 0.021117658*x + -0.002229459*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.011432751 + 0.021566576*x + -0.002105576*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.027393205 + 0.021872871*x + -0.002174845*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.033500217 + 0.012077554*x + -0.001548165*x^2}, size=0.5, xlim=c(0,6), colour='#800C7444') +
  stat_function(fun=function(x){-0.022029305 + 0.015568523*x + -0.001991342*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.053827705 + 0.011180907*x + -0.001431278*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.033191434 + 0.003246487*x + -0.000426906*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.013966838 + 0.004339203*x + -0.0004647*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){-0.031448826 + -0.006390346*x + 0.000483603*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.020152353 + -0.023485325*x + 0.001976154*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.011683515 + -0.032044553*x + 0.002628645*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.04645807 + -0.027919844*x + 0.001982293*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.001717066 + -0.029427476*x + 0.002424643*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.020519863 + -0.008732688*x + 0.001389676*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.016666474 + -0.000610372*x + 0.000657294*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.035720915 + 0.011365366*x + -0.001389324*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){0.001584061 + 0.023960926*x + -0.002109459*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){-0.015995363 + 0.021007355*x + -0.001828975*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.044195593 + 0.021443759*x + -0.0017091*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){-0.036319542 + 0.052472105*x + -0.003171142*x^2}, size=0.5, xlim=c(0,4), colour='#1400E544') +
  stat_function(fun=function(x){-0.012459455 + 0.025667652*x + -0.002161746*x^2}, size=0.5, xlim=c(0,4), colour='#800C7444') +
  stat_function(fun=function(x){-0.042320336 + 0.048739534*x + -0.00292485*x^2}, size=0.5, xlim=c(0,4), colour='#4A06AC44') +
  stat_function(fun=function(x){0.00143156 + -0.003964694*x + 0.0008574*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.046939619 + 0.021907105*x + -0.002626824*x^2}, size=0.5, xlim=c(0,2), colour='#4A06AC44') +
  stat_function(fun=function(x){0.046464844 + 0.016370484*x + -0.002491818*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
  stat_function(fun=function(x){0.022124988 + 0.007186385*x + -0.001828142*x^2}, size=0.5, xlim=c(0,6), colour='#1400E544') +
  stat_function(fun=function(x){0.011055201 + 0.008520216*x + -0.001676888*x^2}, size=0.5, xlim=c(0,6), colour='#4A06AC44') +
  stat_function(fun=function(x){0.048482777 + 0.017583871*x + -0.002622398*x^2}, size=0.5, xlim=c(0,2), colour='#800C7444') +
  stat_function(fun=function(x){0.047044802 + 0.017439179*x + -0.002559794*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.062033316 + 0.025574466*x + -0.003014897*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.001555831 + 0.00461116*x + -0.001580578*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){0.052528883 + 0.019032602*x + -0.002512958*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0.010409216 + 0.012735236*x + -0.001980907*x^2}, size=0.5, xlim=c(0,2), colour='#B6123C44') +
  stat_function(fun=function(x){-0.005363076 + 0.004512983*x + -0.001339912*x^2}, size=0.5, xlim=c(0,2), colour='#EC180444') +
  stat_function(fun=function(x){0.012179776 + 0.012014257*x + -0.001943581*x^2}, size=0.5, xlim=c(0,2), colour='#1400E544') +
#mean lines by plot mani
  #estimated as mean across treatment lines
  #mani1
  stat_function(fun=function(x){(-0.167144 + 0.01842445*x + -0.00178622*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,23), colour='#1400E5') +
  #mani2
  stat_function(fun=function(x){(-0.167144 + 0.01842445*x + -0.00178622*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,22), colour='#4A06AC') +
  #mani3
  stat_function(fun=function(x){(-0.167144 + 0.01842445*x + -0.00178622*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,21), colour='#800C74') +
  #mani4
  stat_function(fun=function(x){(-0.167144 + 0.01842445*x + -0.00178622*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,22), colour='#B6123C') +
  #mani5
  stat_function(fun=function(x){((-0.167144+0.281125) + (0.01842445+0.445834)*x + (-0.00178622-0.0236129)*x^2)*0.1034254 + 0.019179}, size=3, xlim=c(0,22), colour='#EC1804')



evennessInsetPlot <- function() {
  print(evennessPlot)
  
  theme_update(axis.title.x=element_text(size=15, vjust=-0.35, margin=margin(t=5)), axis.text.x=element_text(size=10), axis.title.y=element_text(size=15, angle=90, vjust=0.5, margin=margin(r=5)), axis.text.y=element_text(size=10), plot.title = element_blank())
  
  print(evennessInset, vp=viewport(width = 0.13, height = 0.25, x = 0.80, y = 0.99, just = c("left","top")))
  
  theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=24), axis.title.y=element_text(size=30, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=24))
}

evennessInsetPlot() #export at 1200x1000



###all four variables in one figure
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanInsetPlot(), vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionInsetPlot(), vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessInsetPlot(), vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessInsetPlot(), vp=viewport(layout.pos.row = 2, layout.pos.col = 2))







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



###by resource mani

#mean change
meanResource <- mean%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, water, carbon, yr10)%>%
  gather(key=resource, value=manipulated, nutrients:carbon)%>%
  filter(manipulated!=0)

meanResourcePlot <- ggplot(data=barGraphStats(data=meanResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.5, 0.1), name='Mean Change') +
  coord_cartesian(ylim=c(0, 0.5)) +
  xlab('')+
  annotate('text', x=0.5, y=0.49, label='(a)', size=10, hjust='left')

#dispersion change
dispersionResource <- dispersion%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, water, carbon, yr10)%>%
  gather(key=resource, value=manipulated, nutrients:carbon)%>%
  filter(manipulated!=0)

dispersionResourcePlot <- ggplot(data=barGraphStats(data=dispersionResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.08, 0.1, 0.04), name='Change in Dispersion') +
  coord_cartesian(ylim=c(-0.08, 0.1)) +
  xlab('')+
  annotate('text', x=0.5, y=0.095, label='(b)', size=10, hjust='left')

#richness change
richnessResource <- richness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, water, carbon, yr10)%>%
  gather(key=resource, value=manipulated, nutrients:carbon)%>%
  filter(manipulated!=0)

richnessResourcePlot <- ggplot(data=barGraphStats(data=richnessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(-0.5, 0.2, 0.1), name='Proportion Richness Change') +
  coord_cartesian(ylim=c(-0.22, 0.16)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.145, label='(c)', size=10, hjust='left')

#evenness change
evennessResource <- evenness%>%
  mutate(multi_resource=nutrients+water+carbon)%>%
  filter(multi_resource==1)%>%
  select(site_code, project_name, community_type, treatment, nutrients, water, carbon, yr10)%>%
  gather(key=resource, value=manipulated, nutrients:carbon)%>%
  filter(manipulated!=0)

evennessResourcePlot <- ggplot(data=barGraphStats(data=evennessResource, variable='yr10', byFactorNames=c('resource')), aes(x=resource, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 0.06, 0.01), name='Change in Evenness') +
  coord_cartesian(ylim=c(0, 0.06)) +
  xlab('Resource Manipulated')+
  annotate('text', x=0.5, y=0.058, label='(d)', size=10, hjust='left') 
  
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  

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


#########################################################################################################
#########################################################################################################

#full chains data
chains <- read.csv('fullChains_anpp.csv')%>%
  select(-X)

