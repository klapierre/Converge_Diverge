library(ggplot2)
library(grid)
library(mgcv)
library(lsmeans)
library(codyn)
library(plyr)
library(dplyr)
library(tidyr)

#kim's laptop
setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

#kim's desktop
setwd("C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34, color='black'),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34, color='black'),
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

colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)

##################################################################################
##################################################################################
#experiment information --------------------------------------------------------
expRaw <- read.csv('ExperimentInformation_May2017.csv')

expInfo <- expRaw%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  mutate(irrigation=ifelse(precip>0, 1, 0), drought=ifelse(precip<0, 1, 0))%>%
  summarise(min_year=min(treatment_year), nutrients=mean(nutrients), water=mean(water), carbon=mean(carbon), irrigation=mean(irrigation), drought=mean(drought))

rawData <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\single_resource\\single_resource_data.csv')

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
rawDataAll <- rawData
expInfoSummary <- rawDataAll%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich), anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  summarise(length_mean=mean(experiment_length), length_min=min(experiment_length), length_max=max(experiment_length),
            plot_mani_median=mean(plot_mani), plot_mani_min=min(plot_mani), plot_mani_max=max(plot_mani),
            rrich_mean=mean(rrich), rrich_min=min(rrich), rrich_max=max(rrich),
            anpp_mean=mean(anpp), anpp_min=min(anpp), anpp_max=max(anpp),
            MAP_mean=mean(MAP), MAP_min=min(MAP), MAP_max=max(MAP),
            MAT_mean=mean(MAT), MAT_min=min(MAT), MAT_max=max(MAT))%>%
  gather(variable, estimate)

#treatment info
trtInfo <- rawData%>%
  filter(plot_mani<6, anpp!='NA')%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, resource)%>%
  summarise(experiment_length=mean(experiment_length), plot_mani=mean(plot_mani), rrich=mean(rrich),
            anpp=mean(anpp), MAT=mean(MAT), MAP=mean(MAP))%>%
  ungroup()%>%
  left_join(expInfo)%>%
  mutate(resource_mani=(nutrients+carbon+irrigation+drought), id=1:length(treatment))

################################################################################
################################################################################

# #only run to generate initial chains files
# #raw chains data --------------------------------------------------------
# memory.limit(size=50000)
# chains1 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\single_resource\\single_resource_cholesky_0.csv', comment.char='#')
# chains1 <- chains1[-1:-5000,]
# chains2 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\single_resource\\single_resource_cholesky_1.csv', comment.char='#')
# chains2 <- chains2[-1:-5000,]
# chains3 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\single_resource\\single_resource_cholesky_2.csv', comment.char='#')
# chains3 <- chains3[-1:-5000,]
# chains4 <- read.csv('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\nate_results\\single_resource\\single_resource_cholesky_3.csv', comment.char='#')
# chains4 <- chains4[-1:-5000,]
# 
# chainsCommunity <- rbind(chains1, chains2, chains3, chains4)
# 
# 
# #density plot of chains --------------------------------------------------------
# plot(density(chainsCommunity$mu.1.1))
# plot(density(chainsCommunity$mu.1.2))
# plot(density(chainsCommunity$mu.1.3))


# #get values for overall (mean) lines across levels of plot mani --------------------------------------------------------
# #mean change are the 1's, dispersion are the 2's, richness are the 4's, evenness are the 3's
# chainsCommunity2 <- chainsCommunity%>%
#   select(lp__,
#          #resource intercepts (center digit): 1=drought, 2=irrigation, 3=N, 4=other nutrient, 5=P
#          U.1.1.1, U.2.1.1, U.3.1.1, U.4.1.1,
#          U.1.2.1, U.2.2.1, U.3.2.1, U.4.2.1,
#          U.1.3.1, U.2.3.1, U.3.3.1, U.4.3.1,
#          U.1.4.1, U.2.4.1, U.3.4.1, U.4.4.1,
#          U.1.5.1, U.2.5.1, U.3.5.1, U.4.5.1,
#          #resource linear slopes (center digit): 1=drought, 2=irrigation, 3=N, 4=other nutrient, 5=P
#          U.1.1.2, U.2.1.2, U.3.1.2, U.4.1.2,
#          U.1.2.2, U.2.2.2, U.3.2.2, U.4.2.2,
#          U.1.3.2, U.2.3.2, U.3.3.2, U.4.3.2,
#          U.1.4.2, U.2.4.2, U.3.4.2, U.4.4.2,
#          U.1.5.2, U.2.5.2, U.3.5.2, U.4.5.2,
#          #resource quad slopes (center digit): 1=drought, 2=irrigation, 3=N, 4=other nutrient, 5=P
#          U.1.1.3, U.2.1.3, U.3.1.3, U.4.1.3,
#          U.1.2.3, U.2.2.3, U.3.2.3, U.4.2.3,
#          U.1.3.3, U.2.3.3, U.3.3.3, U.4.3.3,
#          U.1.4.3, U.2.4.3, U.3.4.3, U.4.4.3,
#          U.1.5.3, U.2.5.3, U.3.5.3, U.4.5.3,
#          #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
#          D.1.1.1, D.2.1.1, D.3.1.1, D.4.1.1,
#          D.1.1.2, D.2.1.2, D.3.1.2, D.4.1.2,
#          D.1.1.3, D.2.1.3, D.3.1.3, D.4.1.3,
#          #richness intercept, linear, and quad slopes (center digit): 2=richness
#          D.1.2.1, D.2.2.1, D.3.2.1, D.4.2.1,
#          D.1.2.2, D.2.2.2, D.3.2.2, D.4.2.2,
#          D.1.2.3, D.2.2.3, D.3.2.3, D.4.2.3,
#          #MAP intercept, linear, and quad slopes (center digit): 1=MAP
#          E.1.1.1, E.2.1.1, E.3.1.1, E.4.1.1,
#          E.1.1.2, E.2.1.2, E.3.1.2, E.4.1.2,
#          E.1.1.3, E.2.1.3, E.3.1.3, E.4.1.3,
#          #MAT intercept, linear, and quad slopes (center digit): 2=MAT
#          E.1.2.1, E.2.2.1, E.3.2.1, E.4.2.1,
#          E.1.2.2, E.2.2.2, E.3.2.2, E.4.2.2,
#          E.1.2.3, E.2.2.3, E.3.2.3, E.4.2.3,
#          #overall intercept, linear, and quad slopes (resource=CO2)
#          mu.1.1, mu.2.1, mu.3.1, mu.4.1,
#          mu.1.2, mu.2.2, mu.3.2, mu.4.2,
#          mu.1.3, mu.2.3, mu.3.3, mu.4.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:mu.4.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(lower=median-2*sd, upper=median+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, median=ifelse(diff==-2, 0, median))
# 
# write.csv(chainsCommunity2, 'bayesian_output_summary_single resource_09122017.csv')

chainsCommunity2 <- read.csv('bayesian_output_summary_single resource_09122017.csv')

# #gather the intercepts, linear slopes, and quadratic slopes for all treatments ---------------------------------------------
# #numbers are B.variable.number.parameter (e.g., B.mean.87.slope)
# #variable (second place): 1=mean change, 2=dispersion change, 3=evenness change, 4=richness change
# #parameter (final digit): 1=intercept, 2=linear slope, 3=quad slope
# #set any that are not significant (CI overlaps 0) as 0
# 
# #get mean parameter values across all runs for each experiment, treatment, etc
# chainsFinalMean <- as.data.frame(colMeans(chainsCommunity[,4972:6663]))%>% #may need to delete original four chains dataframes to get this to work
#   add_rownames('parameter')
# names(chainsFinalMean)[names(chainsFinalMean) == 'colMeans(chainsCommunity[, 4972:6663])'] <- 'mean'
# #get sd of parameter values across all runs for each experiment, treatment, etc
# chainsFinalSD <- as.data.frame(colSd(chainsCommunity[,4972:6663]))
# names(chainsFinalSD)[names(chainsFinalSD) == 'colSd(chainsCommunity[, 4972:6663])'] <- 'sd'
# 
# chainsFinal <- cbind(chainsFinalMean, chainsFinalSD)%>%
#   #split names into parts
#   separate(parameter, c('B', 'variable', 'id', 'parameter'))%>%
#   select(-B)%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          id=as.integer(id))%>%
#   #if 95% confidence interval overlaps 0, then set mean to 0
#   mutate(lower=mean-2*sd, upper=mean+2*sd, lower_sign=sign(lower), upper_sign=sign(upper), diff=lower_sign-upper_sign, mean=ifelse(diff==-2, 0, mean))%>%
#   #spread by variable
#   select(variable, id, parameter, mean)%>%
#   spread(key=parameter, value=mean)
# 
# write.csv(chainsFinal, 'bayesian_output_mean sd_single resource_09122017.csv')

chainsFinal <- read.csv('bayesian_output_mean sd_single resource_09122017.csv')

#merge together with experiment list
chainsExperiment <- chainsFinal%>%
  arrange(id)%>%
  left_join(trtInfo, by='id')

chainsEquations <- chainsExperiment%>%
  #get standardized experiment length
  mutate(alt_length=experiment_length - min_year)%>%
  mutate(alt_length=ifelse(alt_length>=9, 8, alt_length))%>%
  mutate(yr9=ifelse(variable=='mean', (intercept+linear*9+quadratic*9^2)*(0.1102327)+(0.2479358),
                    ifelse(variable=='dispersion', (intercept+linear*9+quadratic*9^2)*(0.0764769)+(0.0007942693),
                           ifelse(variable=='evenness', (intercept+linear*9+quadratic*9^2)*(0.09476219)+(0.009980187), (intercept+linear*9+quadratic*9^2)*(0.1781236)+(-0.02631029)))))%>%
  mutate(yr_final=ifelse(variable=='mean', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.1482058)+(0.298937),
                    ifelse(variable=='dispersion', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.08790049)+(-0.000407613),
                           ifelse(variable=='evenness', (intercept+linear*alt_length+quadratic*alt_length^2)*(0.0783639)+(0.01730774), (intercept+linear*alt_length+quadratic*alt_length^2)*(0.2154907)+(-0.0546117)))))%>%
  mutate(curve1='stat_function(fun=function(x){(',
         curve2=' + ',
         curve3='*x + ',
         curve4=ifelse(variable=='mean', '*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,',
                       ifelse(variable=='dispersion', '*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,',
                              ifelse(variable=='evenness', '*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,', '*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,'))),
         curve5='), colour=',
         curve6=') +',
         color=ifelse(resource=='CO2', '#0001E544', ifelse(resource=='drought', '#00A6DD44', ifelse(resource=='irrigation', '#00D56B44', ifelse(resource=='n', '#31CE0044', ifelse(resource=='other_nut', '#C4C60044', '#BF330044'))))),
         curve=paste(curve1, intercept, curve2, linear, curve3, quadratic, curve4, alt_length, curve5, color, curve6, sep='')) 
#need to export this, put quotes around the colors, and copy and paste the curve column back into the ggplot code below
# write.csv(chainsEquations,'plot mani_equations_single resource.csv', row.names=F)



###main figure
# mean change panel --------------------------------------------------------
meanPlot <- ggplot(data=data.frame(x=c(0,0))) + 
  coord_cartesian(ylim=c(0,1))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  ylim(-10,10) +
  xlab('Standardized Year') +
  ylab('Mean Change') +
  annotate('text', x=0, y=1, label='(a)', size=10, hjust='left')

meanPlot <- meanPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0.5260750526682*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.490282746745*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(1.5193225596 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(1.92505883645 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(2.0524468449 + 0.5511237612865*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(1.08157506141 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(1.3248187252 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.559041632028*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0.65988734893 + 0.546203459546695*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.609642048304974 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.6446328655751 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.62961226716476*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.569220580873905*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0.5618057107958*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.562221585406295*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.6812446468992*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.4914432260508*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.631173779905305*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(-0.87803942401105 + 0.695204546225*x + -0.081443104669095*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.4097715966 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(-1.5971575144 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.5190447822 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.4408295979 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.4300662951723*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.6481844725205*x + -0.08131561638321*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.899374652615 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.88205522608 + 0.429758320043755*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.76754929406439 + 0.451022364054922*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.735666701746 + 0.461117626727005*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.79445547148 + 0.4289423495117*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.687524770422225 + 0.60382769178*x + -0.0807649574286975*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.66590283968145 + 0.47680624310218*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.881220733255*x + -0.1158526188352*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.62379057853914 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.73583249648715 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.3604735804 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.14676683578 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.157905016965 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.22518260725 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-1.2643025344 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(1.034876129841 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(1.07985714548 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.6180692877906 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.0015089303465 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.7467585567602 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.643818266561185 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.683971332720885 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(-0.6925229043575 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(-0.53359900683655 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0.5654357795745*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.6473799554309 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0.58340502932461 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0.68352036199995 + 0.42395347158435*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0.76008059612246 + 0.3686109978392*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.69388313254115 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.81983135383273 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.82892495056595 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-1.19882642856 + 0.365145044794555*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.0142652997005 + 0.472154208264*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0.60011687452155*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.9589578236005 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(-0.6311509911687 + 0.715062891985*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0.71900846109605 + 0.4975473494167*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.5662878767188*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0.929418546905 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(1.293359809335 + 0.869301877685*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0.60097308123917 + 0.56996327281805*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0.429262192632215*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.54389428583485 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.797145542206 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.8151110252872 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.74325362494274 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.0935441327 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.011878280175 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(-1.21610223245 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.246527154623 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(-0.80625237876245 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.8261453987885 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,7), colour='#0001E544') +
  stat_function(fun=function(x){(-1.129873086485 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,7), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.6379924944839 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.4138303813 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.134138045547 + 1.07369072645*x + -0.153223627535*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(-1.3930159118 + 0.76376932126945*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.0510010997115 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.8460710688073 + 0.6790520764375*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.046351971617 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.7090442172371*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.626072457770465 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.9952893286415 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.436996292269703*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0.43781578083775*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,6), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1102327)+(0.2479358)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  #last five are the main plot_mani effect lines
  #estimated as mean across treatment lines (plot mani 1-4 staggered by intercept so lines don't overlap)
  #CO2
  stat_function(fun=function(x){(-0.3228245 + 0.209387*x + -0.02306935*x^2)*(0.1102327)+(0.2479358)}, size=3, xlim=c(0,8), colour='#0001E5') +
  #drought
  stat_function(fun=function(x){((-0.3228245-0.10288) + (0.209387+0.05246785)*x + (-0.02306935+0.00516066)*x^2)*(0.1102327)+(0.2479358)}, size=3, xlim=c(0,8), colour='#00A6DD') +
  #irrigation
  stat_function(fun=function(x){((-0.3228245-0.1926415) + (0.209387+0.0229017)*x + (-0.02306935+0.00165397)*x^2)*(0.1102327)+(0.2479358)}, size=3, xlim=c(0,8), colour='#00D56B') +
  #N
  stat_function(fun=function(x){((-0.3228245+0.1005475) + (0.209387+0.04158105)*x + (-0.02306935-0.00110525)*x^2)*(0.1102327)+(0.2479358)}, size=3, xlim=c(0,8), colour='#31CE00') +
  #other_nut
  stat_function(fun=function(x){((-0.3228245-0.02624785) + (0.209387-0.1126835)*x + (-0.02306935-0.1061905)*x^2)*(0.1102327)+(0.2479358)}, size=3, xlim=c(0,8), colour='#C4C600') +
  #P
  stat_function(fun=function(x){((-0.3228245+0.186508) + (0.209387+0.04699845)*x + (-0.02306935-0.008722525)*x^2)*(0.1102327)+(0.2479358)}, size=3, xlim=c(0,8), colour='#BF3300')

# print(meanPlot) #export at 1200x1000


#dispersion panel --------------------------------------------------------
dispersionPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.4))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Dispersion Change') +
  annotate('text', x=0, y=0.4, label='(b)', size=10, hjust='left')

dispersionPlot <- dispersionPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + -0.06964730615337*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0.6927405703161 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(1.596037267769 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.9227595124407 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.9568846915063 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.7669864844659 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.782935560497 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.049771796133 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.8605966572053 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(-2.30323244375 + 1.0483469178495*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-2.857729055 + 1.203243476577*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(-1.7458287092 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.79282511168145 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + -0.079032405268928*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0.8634478304369 + 0.9644503697*x + -0.122083840785*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.746064408245*x + -0.1196196614515*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0.79868373732755 + 0.8758390909*x + -0.123799232815*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + -0.05680935742093*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + -0.055026284489967*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0.8693341177475 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0.955051858103 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(1.1989027059485 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(1.124928251748 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(1.014063246184 + -0.69661354687595*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(1.09060603431391 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.78824022982425 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(1.3421576649285 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(1.7793427467 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(-0.868247246339288 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(-0.938638195832153 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,7), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,7), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0.80740384231135 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.282261167343 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.207294415759 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(1.6270578581415 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.604612047319815*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(-1.5153616266 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,6), colour='#0001E544') +
  stat_function(fun=function(x){(-0.83804467359419 + 0.56959979395395*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.640852194802*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.9825889345555 + 0.5729602091747*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.71844898442145*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0764769)+(0.0007942693)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  #estimated as mean across treatment lines
    #CO2
    stat_function(fun=function(x){(-0.180208 + -0.0182506*x + -0.0077625*x^2)*(0.0764769)+(0.0007942693)}, size=3, xlim=c(0,8), colour='#0001E5') +
    #drought
    stat_function(fun=function(x){((-0.180208+0.115122) + (-0.0182506+0.261739)*x + (-0.0077625-0.03981245)*x^2)*(0.0764769)+(0.0007942693)}, size=3, xlim=c(0,8), colour='#00A6DD') +
    #irrigation
    stat_function(fun=function(x){((-0.180208-0.0163061) + (-0.0182506-0.0274927)*x + (-0.0077625+0.0189937)*x^2)*(0.0764769)+(0.0007942693)}, size=3, xlim=c(0,8), colour='#00D56B') +
    #N
    stat_function(fun=function(x){((-0.180208+0.125067) + (-0.0182506+0.190474)*x + (-0.0077625-0.01741395)*x^2)*(0.0764769)+(0.0007942693)}, size=3, xlim=c(0,8), colour='#31CE00') +
    #other_nut
    stat_function(fun=function(x){((-0.180208+0.268732) + (-0.0182506+0.02267475)*x + (-0.0077625+0.0205392)*x^2)*(0.0764769)+(0.0007942693)}, size=3, xlim=c(0,8), colour='#C4C600') +
    #P
    stat_function(fun=function(x){((-0.180208+0.09573125) + (-0.0182506-0.07421345)*x + (-0.0077625+0.0193059)*x^2)*(0.0764769)+(0.0007942693)}, size=3, xlim=c(0,8), colour='#BF3300')

# print(dispersionPlot) #export at 1200x1000


#richness panel --------------------------------------------------------
richnessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(ylim=c(-1.0,0.8))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.2)) +
  xlab('Standardized Year') +
  ylab('Richness Change') +
  annotate('text', x=0, y=0.8, label='(c)', size=10, hjust='left')

richnessPlot <- richnessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.844623869381 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.912516942338705 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.280481866696 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00A6DD44') +
  stat_function(fun=function(x){(1.1395546955795 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0.9321198615721 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(1.0988838562063 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(1.56300123219 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + -0.39195835279905*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + -0.39032340541669*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0.9533149652065 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + -0.38278330896075*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0.77338384774005 + -0.37539591924504*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + -0.45172069526285*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + -0.43951748775595*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.7728502347644 + -0.5028835827245*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(1.118891537128 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(1.087010982007 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0.455428403155785*x + -0.0756459711482085*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.09789630481865 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(-1.41589207673 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.913878649658 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + -0.048470479874592*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + -0.05134307904227*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + -0.0573201453267*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(1.392856693375 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(1.336152522915 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + -0.42238933650965*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + -0.56343319469745*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + -0.55239829593641*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(1.383305025398 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(2.58431840015 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + -0.69866722467*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + -0.58744085847775*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + -0.5829173774377*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0.962713144715764 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0.870873398899 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,7), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,7), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.7891409766319 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.749960631035785 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(1.261421141776 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(1.27617849043865 + -0.5663734444849*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0.77947331533691 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(1.140850199855 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(-1.102884584878 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.64419610704665*x + -0.082745456717495*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0.57802722655395*x + -0.0783215056356195*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,6), colour='#0001E544') +
  stat_function(fun=function(x){(-1.29764037663 + 0.495927397391955*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.523515488335045*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-2.2149961355 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.3495246342701 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.810059584 + 0*x + 0*x^2)*(0.1781236)+(-0.02631029)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #CO2
  stat_function(fun=function(x){(-0.08033155 + 0.06305445*x + -0.008833285*x^2)*(0.1781236)+(-0.02631029)}, size=3, xlim=c(0,8), colour='#0001E5') +
  #drought
  stat_function(fun=function(x){((-0.08033155-0.2037) + (0.06305445+0.12765)*x + (-0.008833285-0.01021545)*x^2)*(0.1781236)+(-0.02631029)}, size=3, xlim=c(0,8), colour='#00A6DD') +
  #irrigation
  stat_function(fun=function(x){((-0.08033155+0.263938) + (0.06305445-0.13795)*x + (-0.008833285+0.01268175)*x^2)*(0.1781236)+(-0.02631029)}, size=3, xlim=c(0,8), colour='#00D56B') +
  #N
  stat_function(fun=function(x){((-0.08033155+0.198002) + (0.06305445-0.101944)*x + (-0.008833285+0.00613016)*x^2)*(0.1781236)+(-0.02631029)}, size=3, xlim=c(0,8), colour='#31CE00') +
  #other_nut
  stat_function(fun=function(x){((-0.08033155+0.448519) + (0.06305445+0.0856441)*x + (-0.008833285+0.07564105)*x^2)*(0.1781236)+(-0.02631029)}, size=3, xlim=c(0,8), colour='#C4C600') +
  #P
  stat_function(fun=function(x){((-0.08033155+0.134796) + (0.06305445-0.244816)*x + (-0.008833285+0.02778435)*x^2)*(0.1781236)+(-0.02631029)}, size=3, xlim=c(0,8), colour='#BF3300')

# print(richnessPlot) #export at 1200x1000


#evenness panel --------------------------------------------------------
evennessPlot <- ggplot(data=data.frame(x=c(0,0))) +
  # coord_cartesian(ylim=c(-0.05,0.35))  +
  scale_x_continuous(limits=c(0,8), breaks=seq(0,8,1), labels=seq(1,9,1)) +
  # scale_y_continuous(limits=c(-2,2), breaks=seq(-2,2,0.1)) +
  xlab('Standardized Year') +
  ylab('Evenness Change') +
  annotate('text', x=0, y=0.6, label='(d)', size=10, hjust='left')

evennessPlot <- evennessPlot + 
  #below are the individual treatment lines
  stat_function(fun=function(x){(1.2141325770105 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(1.197006149564 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + -0.088033390657375*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + -0.085589761539055*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(-1.49875748491 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(-1.4054960729 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.3277971218899 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0.98391423166055 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(1.02278679107962 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#C4C60044') +
  stat_function(fun=function(x){(-0.87680749682955 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.5700141941928*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.67600644482665*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0.4916909845727*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.434125889095525*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.46840001511432*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.454956928540111*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.418741851508*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.51769500003052*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.4998247857985*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.847378204830565 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.8905633202275 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.9174059490234 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(-1.218800204795 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-1.2292212225 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0.984207114962 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(-0.63950266029245 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.6548040488432 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.9776406449609 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(-0.996892015109 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.5929850579835*x + -0.093750738143*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(-3.261400496 + 1.19469078211*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(-3.3172226775 + 1.1734342726395*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.81522398974195 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(-0.82106445834775 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 1.014960541825*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0.96792046952*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.763049674369*x + -0.11072328559135*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0.60953089980705*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,7), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,7), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-0.82756293142868 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,5), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,3), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0.8893221567725*x + -0.1245520237325*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,8), colour='#0001E544') +
  stat_function(fun=function(x){(-0.79071863132842 + 0.82509202548875*x + -0.11672504134142*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + -0.07745697769658*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(-1.3071336036155 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0.9645984015848 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00A6DD44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#00D56B44') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.90212601657375 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#31CE0044') +
  stat_function(fun=function(x){(-0.8309756702585 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(-0.78903002068005 + 0.47824409857107*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,4), colour='#BF330044') +
  stat_function(fun=function(x){(2.4201899999 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#BF330044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,6), colour='#0001E544') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  stat_function(fun=function(x){(0 + 0*x + 0*x^2)*(0.0783639)+(0.009980187)}, size=0.5, xlim=c(0,2), colour='#31CE0044') +
  #mean lines by plot mani
  #estimated as mean across treatment lines
  #CO2
  stat_function(fun=function(x){(-0.216373 + -0.0258594*x + -0.0037162*x^2)*(0.0783639)+(0.009980187)}, size=3, xlim=c(0,8), colour='#0001E5') +
  #drought
  stat_function(fun=function(x){((-0.216373+0.062114) + (-0.0258594+0.03989015)*x + (-0.0037162-0.02475195)*x^2)*(0.0783639)+(0.009980187)}, size=3, xlim=c(0,8), colour='#00A6DD') +
  #irrigation
  stat_function(fun=function(x){((-0.216373-0.2994262) + (-0.0258594+0.1175841)*x + (-0.0037162-0.01336313)*x^2)*(0.0783639)+(0.009980187)}, size=3, xlim=c(0,8), colour='#00D56B') +
  #N
  stat_function(fun=function(x){((-0.216373-0.077241) + (-0.0258594+0.0986766)*x + (-0.0037162-0.0155275)*x^2)*(0.0783639)+(0.009980187)}, size=3, xlim=c(0,8), colour='#31CE00') +
  #other_nut
  stat_function(fun=function(x){((-0.216373+0.620346) + (-0.0258594-0.2097594)*x + (-0.0037162-0.1841572)*x^2)*(0.0783639)+(0.009980187)}, size=3, xlim=c(0,8), colour='#C4C600') +
  #P
  stat_function(fun=function(x){((-0.216373-0.061326) + (-0.0258594+0.2222626)*x + (-0.0037162-0.0412204)*x^2)*(0.0783639)+(0.009980187)}, size=3, xlim=c(0,8), colour='#BF3300')

# print(evennessPlot) #export at 1200x1000


#print all plots together --------------------------------------------------------
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersionPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(richnessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(evennessPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))
#export at 2400 x 2000


# ##summary stats from bayesian output --------------------------------------------------------
# #gather summary stats needed and relabel them
# chainsCommunitySummary <- chainsCommunity%>%
#   select(#resource_type intercepts (center digit): 1=drought, 2=irr, 3=N, 4=othernut, 5=p
#     U.1.1.1, U.2.1.1, U.3.1.1, U.4.1.1,
#     U.1.2.1, U.2.2.1, U.3.2.1, U.4.2.1,
#     U.1.3.1, U.2.3.1, U.3.3.1, U.4.3.1,
#     U.1.4.1, U.2.4.1, U.3.4.1, U.4.4.1,
#     U.1.5.1, U.2.5.1, U.3.5.1, U.4.5.1,
#     #plot_mani linear slopes (center digit): 1=drought, 2=irr, 3=N, 4=othernut, 5=p
#     U.1.1.2, U.2.1.2, U.3.1.2, U.4.1.2,
#     U.1.2.2, U.2.2.2, U.3.2.2, U.4.2.2,
#     U.1.3.2, U.2.3.2, U.3.3.2, U.4.3.2,
#     U.1.4.2, U.2.4.2, U.3.4.2, U.4.4.2,
#     U.1.5.2, U.2.5.2, U.3.5.2, U.4.5.2,
#     #plot_mani quad slopes (center digit): 1=drought, 2=irr, 3=N, 4=othernut, 5=p
#     U.1.1.3, U.2.1.3, U.3.1.3, U.4.1.3,
#     U.1.2.3, U.2.2.3, U.3.2.3, U.4.2.3,
#     U.1.3.3, U.2.3.3, U.3.3.3, U.4.3.3,
#     U.1.4.3, U.2.4.3, U.3.4.3, U.4.4.3,
#     U.1.5.3, U.2.5.3, U.3.5.3, U.4.5.3,
#     #ANPP intercept, linear, and quad slopes (center digit): 1=anpp
#     D.1.1.1, D.2.1.1, D.3.1.1, D.4.1.1,
#     D.1.1.2, D.2.1.2, D.3.1.2, D.4.1.2,
#     D.1.1.3, D.2.1.3, D.3.1.3, D.4.1.3,
#     #richness intercept, linear, and quad slopes (center digit): 2=richness
#     D.1.2.1, D.2.2.1, D.3.2.1, D.4.2.1,
#     D.1.2.2, D.2.2.2, D.3.2.2, D.4.2.2,
#     D.1.2.3, D.2.2.3, D.3.2.3, D.4.2.3,
#     #MAP intercept, linear, and quad slopes (center digit): 1=MAP
#     E.1.1.1, E.2.1.1, E.3.1.1, E.4.1.1,
#     E.1.1.2, E.2.1.2, E.3.1.2, E.4.1.2,
#     E.1.1.3, E.2.1.3, E.3.1.3, E.4.1.3,
#     #MAT intercept, linear, and quad slopes (center digit): 2=MAT
#     E.1.2.1, E.2.2.1, E.3.2.1, E.4.2.1,
#     E.1.2.2, E.2.2.2, E.3.2.2, E.4.2.2,
#     E.1.2.3, E.2.2.3, E.3.2.3, E.4.2.3,
#     #overall intercept, linear, and quad slopes
#     mu.1.1, mu.2.1, mu.3.1, mu.4.1,
#     mu.1.2, mu.2.2, mu.3.2, mu.4.2,
#     mu.1.3, mu.2.3, mu.3.3, mu.4.3)%>%
#   gather(key=parameter, value=value, U.1.1.1:mu.4.3)%>%
#   group_by(parameter)%>%
#   summarise(median=median(value), sd=sd(value))%>%
#   mutate(CI=sd*2)%>%
#   separate(parameter, c('level', 'variable', 'predictor', 'parameter'))%>%
#   mutate(parameter=ifelse(level=='mu', predictor, parameter), predictor=ifelse(level=='mu', 'overall', predictor))%>%
#   #rename parts to be more clear
#   mutate(variable=ifelse(variable==1, 'mean', ifelse(variable==2, 'dispersion', ifelse(variable==3, 'evenness', 'richness'))),
#          parameter=ifelse(parameter==1, 'intercept', ifelse(parameter==2, 'linear', 'quadratic')),
#          predictor=ifelse(level=='D'&predictor==1, 'ANPP', ifelse(level=='D'&predictor==2, 'rrich', ifelse(level=='E'&predictor==1, 'MAP', ifelse(level=='E'&predictor==2, 'MAT', ifelse(level=='U'&predictor==1, 'drought', ifelse(level=='U'&predictor==2, 'irr', ifelse(level=='U'&predictor==3, 'N', ifelse(level=='U'&predictor==4, 'other nut', ifelse(level=='U'&predictor==5, 'P', 'overall'))))))))))%>%
#   select(level, parameter, variable, predictor, predictor, median, sd, CI)
# 
# write.csv(chainsCommunitySummary, 'bayesian_output_summary_final plots_single resource_09122017.csv')
chainsCommunitySummary <- read.csv('bayesian_output_summary_final plots_single resource_09122017.csv')

chainsCommunityOverall <- chainsCommunitySummary%>%
  filter(predictor=='overall')%>%
  mutate(overall=median)%>%
  select(parameter, variable, overall)%>%
  left_join(chainsCommunitySummary, by=c('parameter', 'variable'))%>%
  mutate(overall=ifelse(predictor=='overall', 0, overall))%>%
  mutate(median_corrected=median+overall)




###overall responses from bayesian output --------------------------------------------------------
meanOverallPlot <- ggplot(data=subset(chainsCommunityOverall, variable=='mean' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.5, 0.5, 0.5)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('\nMean Change') +
  annotate('text', x=3.45, y=-0.8, label='(a)', size=10, hjust='left')

dispersionOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='dispersion' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.18), breaks=seq(-0.2, 0.2, 0.2)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Dispersion\nChange') +
  annotate('text', x=3.45, y=-0.3, label='(b)', size=10, hjust='left')

richnessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='richness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.3, 0.5), breaks=seq(-0.3, 0.5, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Richness Change') +
  annotate('text', x=3.45, y=-0.3, label='(c)', size=10, hjust='left')

evennessOverallPlot <- ggplot(data=subset(chainsCommunitySummary, variable=='evenness' & predictor=='overall'), aes(x=parameter, y=median)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_y_continuous(limits=c(-0.4, 0.25), breaks=seq(-0.3, 0.3, 0.3)) +
  scale_x_discrete(limits=c('quadratic', 'linear', 'intercept'),
                   labels=c('Quadratic Slope', 'Linear Slope', 'Intercept')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=28, vjust=2, margin=margin(b=15))) +
  geom_hline(aes(yintercept=0)) +
  coord_flip() +
  ggtitle('Evenness\nChange') +
  annotate('text', x=3.45, y=-0.4, label='(d)', size=10, hjust='left')

pushViewport(viewport(layout=grid.layout(1,4)))
print(evennessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 4))
print(richnessOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(dispersionOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(meanOverallPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
#export at 2400x500



#mean plots --------------------------------------------------------
meanIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.2)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'P', 'other_nut', 'N', 'irrigation', 'drought'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', 'P', 'other nut', 'N', 'irr', 'drought')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  # ylim(-1.15, 1.15) +
  geom_hline(aes(yintercept=0)) +
  coord_flip()

meanSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

meanQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='mean'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#dispersion plots --------------------------------------------------------
dispersionIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

dispersionSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

dispersionQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='dispersion'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#richness plots --------------------------------------------------------
richnessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

richnessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

richnessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='richness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#evenness plots --------------------------------------------------------
evennessIntPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='intercept'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.15, 1.15) +
  coord_flip()

evennessSlopePlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='linear'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  ylim(-1.2, 1.2) +
  coord_flip()

evennessQuadPlot <- ggplot(data=subset(chainsCommunitySummary, parameter=='quadratic'&variable=='evenness'&predictor!='overall'), aes(x=predictor, y=median)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=median-CI, ymax=median+CI, width=0.1)) +
  scale_x_discrete(limits=c('MAT', 'MAP', 'ANPP', 'rrich', 'plot mani 5', 'plot mani 4', 'plot mani 3', 'plot mani 2'),
                   labels=c('MAT', 'MAP', 'ANPP', 'Gamma Diversity', '5 Manipulations', '4 Manipulations', '3 Manipulations', '2 Manipulations')) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_vline(aes(xintercept=4.5), linetype='dashed') +
  geom_vline(aes(xintercept=2.5), linetype='dashed') +
  geom_hline(aes(yintercept=0)) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, 0.1), limits=c(-0.13,0.1)) +
  coord_flip()

#plot all together --------------------------------------------------------
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





# ###by resource mani - model output--------------------------------------------------------
# #still need to calculate the proportion of chains where x resource response was greater than y resource response
# trtDetail <- expRaw%>%
#   select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
#   group_by(site_code, project_name, community_type, treatment)%>%
#   summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
#   mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))
# 
# chainsTrt <- chainsEquations%>%
#   select(variable, site_code, project_name, community_type, treatment, intercept, linear, quadratic, yr9, yr_final, plot_mani, rrich, anpp, MAT, MAP, min_year, experiment_length, alt_length)%>%
#   left_join(trtDetail)
# 
# resourceMani <- chainsTrt%>%
#   mutate(resource_mani=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
#                               ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
#                                      ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
#                                             ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
#                                                    ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'nuts:CO2',
#                                                           ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'nuts:dro',
#                                                                  ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'nuts:irr',
#                                                                         ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'CO2:dro',
#                                                                                ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'CO2:irr',
#                                                                                       ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'nuts:CO2:dro', 
#                                                                                              ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr', 'other'))))))))))))
# 
# #plot model output by resource manipulated at final year of each experiment (varies by experiment) ---------------------------
# meanResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='mean'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.20), name='Mean Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(0, 0.6), xlim=c(1,4)) +
#   xlab('')+
#   annotate('text', x=0.5, y=0.60, label='(a)', size=12, hjust='left') +
#   annotate('text', x=1, y=0.42, label='a*', size=10) +
#   annotate('text', x=2, y=0.47, label='ab*', size=10) +
#   annotate('text', x=3, y=0.35, label='b*', size=10) +
#   annotate('text', x=4, y=0.35, label='b*', size=10)
# 
# dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='dispersion'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.01), name='Dispersion Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.05, 0.03), xlim=c(1,4)) +
#   xlab('') +
#   annotate('text', x=0.5, y=0.03, label='(b)', size=12, hjust='left')
# 
# richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='richness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Richness Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.15, 0.06), xlim=c(1,4)) +
#   xlab('') +
#   annotate('text', x=0.5, y=0.06, label='(c)', size=12, hjust='left') +
#   annotate('text', x=1, y=-0.13, label='ab*', size=10) +
#   annotate('text', x=2, y=-0.10, label='a*', size=10) +
#   annotate('text', x=3, y=-0.05, label='b', size=10) +
#   annotate('text', x=4, y=-0.12, label='ab*', size=10)
# 
# evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=subset(resourceMani, variable=='evenness'&resource_mani!='other'&resource_mani!='nuts:CO2'&resource_mani!='nuts:dro'&resource_mani!='nuts:irr'&resource_mani!='CO2:dro'&resource_mani!='CO2:irr'&resource_mani!='nuts:CO2:dro'&resource_mani!='nuts:CO2:irr'), variable='yr_final', byFactorNames=c('resource_mani')), aes(x=resource_mani, y=mean)) +
#   geom_bar(stat="identity", fill='white', color='black') +
#   geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
#   scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Evenness Change') +
#   scale_x_discrete(limits=c('nuts', 'CO2', 'irrigation', 'drought'),
#                    labels=c('+nutrients', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O)) +
#   coord_cartesian(ylim=c(-0.03, 0.07), xlim=c(1,4)) +
#   xlab('') +
#   annotate('text', x=0.5, y=0.07, label='(d)', size=12, hjust='left') +
#   annotate('text', x=1, y=0.055, label='a*', size=10) +
#   annotate('text', x=3, y=0.05, label='ab*', size=10) +
#   annotate('text', x=4, y=0.02, label='b', size=10)
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600
# 
# 
# ###by magnitude of resource manipulated---------------------------------
# #N addition
# meanNPlotFinal <- ggplot(data=subset(resourceMani, variable=='mean'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=1, label='(a)', size=12, hjust='left')
# 
# dispersionNPlotFinal <- ggplot(data=subset(resourceMani, variable=='dispersion'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=0.5, label='(b)', size=12, hjust='left')
# 
# richnessNPlotFinal <- ggplot(data=subset(resourceMani, variable=='richness'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=2, label='(c)', size=12, hjust='left')
# 
# evennessNPlotFinal <- ggplot(data=subset(resourceMani, variable=='evenness'&n>0), aes(x=n, y=yr_final)) +
#   geom_point(size=5) +
#   coord_trans(x="log") +
#   scale_y_continuous(breaks=seq(-0.35, 0.5, 0.25), name='Evenness Change') +
#   xlab('N added (gm-2)') +
#   annotate('text', x=0.1, y=0.5, label='(d)', size=12, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600
# 
# #H2O change
# meanPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='mean'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=1, label='(a)', size=12, hjust='left')
# 
# dispersionPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='dispersion'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=0.5, label='(b)', size=12, hjust='left')
# 
# richnessPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='richness'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=2, label='(c)', size=12, hjust='left')
# 
# evennessPrecipPlotFinal <- ggplot(data=subset(resourceMani, variable=='evenness'&precip!=0), aes(x=precip, y=yr_final)) +
#   geom_point(size=5) +
#   # scale_y_continuous(breaks=seq(-0.35, 0.5, 0.25), name='Evenness Change') +
#   xlab('H2O added/removed (mm)') +
#   annotate('text', x=0.1, y=0.5, label='(d)', size=12, hjust='left')
# 
# pushViewport(viewport(layout=grid.layout(2,2)))
# print(meanPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
# print(dispersionPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
# print(richnessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(evennessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
# #export at 1800 x 1600


###by resource mani - raw data--------------------------------------------------------
trtDetail <- expRaw%>%
  select(site_code, project_name, community_type, treatment, n, p, k, CO2, precip)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(n=mean(n), p=mean(p), k=mean(k), CO2=mean(CO2), precip=mean(precip))%>%
  mutate(drought=ifelse(precip<0, precip, 0), irrigation=ifelse(precip>0, precip, 0))

rawTrt <- rawData%>%
  filter(treatment_year==9|treatment_year==experiment_length)%>%
  select(site_code, project_name, community_type, treatment, plot_mani, rrich, anpp, MAT, MAP, experiment_length, treatment_year, mean_change, dispersion_change, SimpEven_change, S_PC)%>%
  left_join(trtDetail)%>%
  mutate(resource_mani=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
                              ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
                                     ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
                                            ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
                                                   ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'nuts:CO2',
                                                          ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'nuts:dro',
                                                                 ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'nuts:irr',
                                                                        ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'CO2:dro',
                                                                               ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'CO2:irr',
                                                                                      ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'nuts:CO2:dro', 
                                                                                             ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'nuts:CO2:irr','other'))))))))))))%>%
  mutate(resource_mani_combo=ifelse((n+p+k)>0&CO2==0&drought==0&irrigation==0, 'nuts',
                                    ifelse((n+p+k)==0&CO2>0&drought==0&irrigation==0, 'CO2',
                                           ifelse((n+p+k)==0&CO2==0&drought<0&irrigation==0, 'drought',
                                                  ifelse((n+p+k)==0&CO2==0&drought==0&irrigation>0, 'irrigation',
                                                         ifelse((n+p+k)>0&CO2>0&drought==0&irrigation==0, 'multiple',
                                                                ifelse((n+p+k)>0&CO2==0&drought<0&irrigation==0, 'multiple',
                                                                       ifelse((n+p+k)>0&CO2==0&drought==0&irrigation>0, 'multiple',
                                                                              ifelse((n+p+k)==0&CO2>0&drought<0&irrigation==0, 'multiple',
                                                                                     ifelse((n+p+k)==0&CO2>0&drought==0&irrigation>0, 'multiple',
                                                                                            ifelse((n+p+k)>0&CO2>0&drought<0&irrigation==0,'multiple', 
                                                                                                   ifelse((n+p+k)>0&CO2>0&drought==0&irrigation>0,'multiple','other'))))))))))))

#plot raw data by resource manipulated at final year of each experiment (varies by experiment) ---------------------------
meanResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='mean_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.10), name='Mean Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(0, 0.45), xlim=c(1,6)) +
  xlab('')+
  annotate('text', x=0.5, y=0.45, label='(a)', size=12, hjust='left') +
  annotate('text', x=1, y=0.35, label='b*', size=10) +
  annotate('text', x=2, y=0.42, label='a*', size=10) +
  annotate('text', x=3, y=0.335, label='b*', size=10) +
  annotate('text', x=4, y=0.35, label='b*', size=10) +
  annotate('text', x=5, y=0.355, label='b*', size=10) +
  annotate('text', x=6, y=0.33, label='b*', size=10)

dispersionResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='dispersion_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.02), name='Dispersion Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.07, 0.03), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.03, label='(b)', size=12, hjust='left')

richnessResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='S_PC', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.1), name='Richness Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.2, 0.15), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.15, label='(c)', size=12, hjust='left') +
  annotate('text', x=1, y=0.07, label='a', size=10) +
  annotate('text', x=2, y=-0.185, label='b*', size=10) +
  annotate('text', x=3, y=0.14, label='ab', size=10) +
  annotate('text', x=4, y=0.125, label='a*', size=10) +
  annotate('text', x=5, y=-0.115, label='ab', size=10) +
  annotate('text', x=6, y=0.09, label='a', size=10)

evennessResourcePlotFinal <- ggplot(data=barGraphStats(data=rawTrt, variable='SimpEven_change', byFactorNames=c('resource_mani_combo')), aes(x=resource_mani_combo, y=mean)) +
  geom_bar(stat="identity", fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(-5, 5, 0.05), name='Evenness Change') +
  scale_x_discrete(limits=c('other', 'nuts', 'CO2', 'irrigation', 'drought', 'multiple'),
                   labels=c('non-', '+nuts', '+' ~CO[2], '+' ~H[2]*O, '-' ~H[2]*O, 'multiple')) +
  coord_cartesian(ylim=c(-0.09, 0.14), xlim=c(1,6)) +
  xlab('') +
  annotate('text', x=0.5, y=0.14, label='(d)', size=12, hjust='left') +
  annotate('text', x=1, y=0.035, label='ab', size=10) +
  annotate('text', x=2, y=0.058, label='a*', size=10) +
  annotate('text', x=3, y=0.128, label='ab', size=10) +
  annotate('text', x=4, y=-0.05, label='b', size=10) +
  annotate('text', x=5, y=-0.055, label='ab', size=10) +
  annotate('text', x=6, y=-0.055, label='b*', size=10)

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionResourcePlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessResourcePlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 2400 x 1600


###by magnitude of resource manipulated---------------------------------
#N addition
meanNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=mean_change)) +
  geom_point(size=5) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=0.4, y=1, label='(a)', size=12, hjust='left')

dispersionNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=dispersion_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-0.4, 0.4, 0.2), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=0.4, y=0.4, label='(b)', size=12, hjust='left')

richnessNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=S_PC)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-1,2,0.5), name='Richness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=1.5, label='(c)', size=12, hjust='left')

evennessNPlotFinal <- ggplot(data=subset(rawTrt, n>0), aes(x=n, y=SimpEven_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_x_log10() +
  scale_y_continuous(breaks=seq(-0.2, 0.6, 0.2), name='Evenness Change') +
  xlab(expression(paste('N added (g', m^-2, ')'))) +
  annotate('text', x=0.4, y=0.6, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionNPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessNPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600

#H2O change
meanPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=mean_change)) +
  geom_point(size=5) +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=-80, y=0.65, label='(a)', size=12, hjust='left')

dispersionPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=dispersion_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=-80, y=0.2, label='(b)', size=12, hjust='left')

richnessPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=S_PC)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks=seq(-1,2,0.25), name='Richness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.6, label='(c)', size=12, hjust='left')

evennessPrecipPlotFinal <- ggplot(data=subset(rawTrt, precip!=0), aes(x=precip, y=SimpEven_change)) +
  geom_point(size=5) +
  geom_hline(yintercept=0) +
  scale_y_continuous(breaks=seq(-0.4, 0.5, 0.1), name='Evenness Change') +
  xlab(expression(paste(H[2], 'O deviation from ambient (%)'))) +
  annotate('text', x=-80, y=0.2, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionPrecipPlotFinal, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessPrecipPlotFinal, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600


###for factorial experiments, are responses (agnostic to factor identity) additive?
#raw data--------------------------------------------------------
rawFactor <- expRaw%>%
  select(site_code, project_name, community_type, treatment, factorial)%>%
  filter(factorial==1)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  unique()%>%
  left_join(rawData)%>%
  filter(!is.na(mean_change))%>%
  filter(treatment_year==9|treatment_year==experiment_length)%>%
  filter(plot_mani<6)

#plot factorial responses for each metric
meanFactorial <- ggplot(barGraphStats(data=rawFactor, variable="mean_change", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
  geom_bar(stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
  scale_y_continuous(breaks=seq(0, 1, 0.20), name='Mean Change') +
  xlab('') +
  annotate('text', x=0.5, y=0.65, label='(a)', size=12, hjust='left')

dispersionFactorial <- ggplot(barGraphStats(data=rawFactor, variable="dispersion_change", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
  geom_bar(stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
  scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1), name='Dispersion Change') +
  xlab('') +
  annotate('text', x=0.5, y=0.38, label='(b)', size=12, hjust='left')

richnessFactorial <- ggplot(barGraphStats(data=rawFactor, variable="S_PC", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
  geom_bar(stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
  scale_y_continuous(breaks=seq(-1,2,0.25), name='Richness Change') +
  xlab('Number of Factors') +
  annotate('text', x=0.5, y=0.6, label='(c)', size=12, hjust='left')

evennessFactorial <- ggplot(barGraphStats(data=rawFactor, variable="SimpEven_change", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
  geom_bar(stat="identity", color='black', fill='white') +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=0.1) +
  scale_y_continuous(breaks=seq(-0.4, 0.5, 0.1), name='Evenness Change') +
  xlab('Number of Factors') +
  annotate('text', x=0.5, y=0.3, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanFactorial, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionFactorial, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(evennessFactorial, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(richnessFactorial, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
#export at 1800 x 1600




#look at number replicates for dispersion results (all factors actually) -- doesn't make a difference --------------------------------------
reps <- read.csv('SpeciesRelativeAbundance_Dec2016.csv')%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year, plot_id)%>%
  summarise(mean=mean(relcov))%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment, calendar_year)%>%
  summarise(rep_num=n())%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarise(rep_num=mean(rep_num))

dispersionReps <- rawTrt%>%
  left_join(reps, by=c('site_code', 'project_name', 'community_type', 'treatment'))

meanRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=mean_change)) +
  geom_point() +
  xlab('Number of Relicates') +
  ylab('Mean Change') +
  scale_x_continuous(breaks=seq(0,50,10)) +
  coord_cartesian(xlim=c(2,45)) +
  annotate('text', x=0, y=1, label='(a)', size=12, hjust='left')
dispersionRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=dispersion_change)) +
  geom_point() +
  geom_hline(yintercept=0) +
  xlab('Number of Relicates') +
  ylab('Dispersion Change') +
  scale_x_continuous(breaks=seq(0,50,10)) +
  coord_cartesian(xlim=c(2,45)) +
  annotate('text', x=0, y=0.35, label='(b)', size=12, hjust='left')
richnessRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=S_PC)) +
  geom_point() +
  geom_hline(yintercept=0) +
  xlab('Number of Relicates') +
  ylab('Richness Change') +
  scale_x_continuous(breaks=seq(0,50,10)) +
  coord_cartesian(xlim=c(2,45)) +
  annotate('text', x=0, y=1.4, label='(c)', size=12, hjust='left')
evennessRepPlot <- ggplot(data=dispersionReps, aes(x=rep_num, y=SimpEven_change)) +
  geom_point() +
  geom_hline(yintercept=0) +
  xlab('Number of Relicates') +
  ylab('Evenness Change') +
  scale_x_continuous(breaks=seq(0,50,10)) +
  coord_cartesian(xlim=c(2,45)) +
  annotate('text', x=0, y=0.6, label='(d)', size=12, hjust='left')

pushViewport(viewport(layout=grid.layout(2,2)))
print(meanRepPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(dispersionRepPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(richnessRepPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(evennessRepPlot, vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
#export at 1800 x 1600







###look at five factor manipulations for mean change --------------------------------------------------------
#compare any four factor without N to five factor with N - using model data final year--------------------------
expRawMean <- expRaw%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarise(n=mean(n), herb_removal=mean(herb_removal), plant_mani=mean(plant_mani))

meanCompare <- subset(chainsEquations, project_name=='e001'|project_name=='e002'|site_code=='NIN'|site_code=='TRA')%>%
  select(variable, site_code, project_name, community_type, treatment, plot_mani, yr9)%>%
  left_join(expRawMean, by=c('site_code', 'project_name', 'community_type', 'treatment', 'plot_mani'), all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))

#plot without N at four factors, with N at five factors
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&variable=='mean'), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none')
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&herb_removal>0&variable=='mean'), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none')
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plant_mani>0&variable=='mean'), variable='yr9', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
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

#compare any four factor without N to five factor with N - using raw data final year--------------------------
expRawMean <- expRaw%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarise(n=mean(n), herb_removal=mean(herb_removal), plant_mani=mean(plant_mani))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, plot_mani, n, herb_removal, plant_mani)

meanCompare <- subset(rawTrt, project_name=='e001'|project_name=='e002'|site_code=='NIN'|site_code=='TRA')%>%
  left_join(expRawMean, all=F)%>%
  mutate(n_mani=ifelse(n>0, 1, 0))


#plot without N at four factors, with N at five factors--------------------------
compareNPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6), variable='mean_change', byFactorNames=c('plot_mani', 'n_mani')), aes(x=interaction(plot_mani, n_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='Mean Change') +
  scale_x_discrete(labels=c('4 factor\n-N', '4 factor\n+N', '5 factor\n+N')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(a) Nitrogen Comparison', size=10, hjust='left') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.4, label='a*', size=10) +
  annotate('text', x=2, y=0.57, label='a*', size=10) +
  annotate('text', x=3, y=0.8, label='b*', size=10)
compareHerbPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6), variable='mean_change', byFactorNames=c('plot_mani', 'herb_removal')), aes(x=interaction(plot_mani, herb_removal), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(limits=c('4.0', '4.1', '5.1'), labels=c('4 factor\n-excl.', '4 factor\n+excl.', '5 factor\n+excl.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('Number of Factors Manipulated') +
  annotate('text', x=0.5, y=1, label='(b) Herbivore Removal Comparison', size=10, hjust='left') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.39, label='*', size=10) +
  annotate('text', x=2, y=0.6, label='*', size=10) +
  annotate('text', x=3, y=0.59, label='*', size=10)
comparePlantPlot <- ggplot(data=barGraphStats(data=subset(meanCompare, plot_mani>3&plot_mani<6), variable='mean_change', byFactorNames=c('plot_mani', 'plant_mani')), aes(x=interaction(plot_mani, plant_mani), y=mean, fill=as.factor(plot_mani))) +
  geom_bar(stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-1.96*se, ymax=mean+1.96*se, width=0.2)) +
  scale_y_continuous(breaks=seq(0, 1.0, 0.2), name='') +
  scale_x_discrete(limits=c('4.0', '4.1', '5.1'), labels=c('4 factor\n-manip', '4 factor\n+manip.', '5 factor\n+manip.')) +
  coord_cartesian(ylim=c(0,1)) +
  scale_fill_manual(values=c('white', 'grey')) +
  xlab('') +
  annotate('text', x=0.5, y=1, label='(c) Plant Manipulation Comparison', size=10, hjust='left') +
  theme(legend.position='none') +
  annotate('text', x=1, y=0.39, label='*', size=10) +
  annotate('text', x=2, y=0.61, label='*', size=10) +
  annotate('text', x=3, y=0.69, label='*', size=10)

pushViewport(viewport(layout=grid.layout(1,3)))
print(compareNPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(comparePlantPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 3))
print(compareHerbPlot, vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
#export at 2400x1200