library(ggplot2)
library(plyr)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

setwd('C:\\Users\\Kim\\Dropbox\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

#read in Meghan's data info
info <- read.csv('exp_info072015.csv')

#read in full change in mean and dispersion dataset
full <- read.csv('dispersion_and_means_press_experiments_with_exp_info_03232015.csv')

#make dataset indicating what control treatments are, plus all other relevant information
controls <- full[,-c(1:5, 7:8, 10:74)]
controlsAggregate <- aggregate(plot_mani.x ~ label, mean, data=full)
exptLength <- aggregate(trt.year ~ label, max, data=full)
names(exptLength)[names(exptLength)=='trt.year'] <- 'expt_length'
controlsInfo <- merge(controlsAggregate, exptLength)

#using full dataset, get slopes of means through time for each treatment independantly
slopesMean <- ddply(subset(full, mean.disp=='mean'), c('expt', 'label'), function(x) {
  model <- lm(dist~trt.year, data=x)
  coef(model)
})
names(slopesMean)[names(slopesMean)=="trt.year"] <- "slope"
  
#merge slopes data with controls and extra info
slopesMeanControls <- merge(slopesMean, controlsInfo)
slopesMeanInfo <- merge(slopesMeanControls, info)

# summary(slopeMeansModel <- lm(slope ~ plot_mani.x, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=plot_mani.x, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ expt_length, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=expt_length, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ MAP, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=MAP, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ ANPP, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=ANPP, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ species_num, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=species_num, y=slope)) +
#   geom_point()

#using full dataset, get slopes of dispersion through time for each treatment independantly
slopesDispersion <- ddply(subset(full, mean.disp=='disp'), c('expt', 'label'), function(x) {
  model <- lm(dist~trt.year, data=x)
  coef(model)
})
names(slopesDispersion)[names(slopesDispersion)=="trt.year"] <- "slope"

#merge slopes data with controls and extra info
slopesDispersionInfo <- merge(slopesDispersion, controlsAggregate)

#get treatment and control effects as columns
slopesDispersionInfoTrt <- subset(slopesDispersionInfo, subset=(plot_mani.x!=0))
slopesDispersionInfoCtl <- subset(slopesDispersionInfo, subset=(plot_mani.x==0))
slopesDispersionInfoCtlSlim <- slopesDispersionInfoCtl[,-c(1, 3, 5)]
names(slopesDispersionInfoCtlSlim)[names(slopesDispersionInfoCtlSlim)=='slope'] <- 'ctl_slope'
slopesDispersionInfoRR <- merge(slopesDispersionInfoTrt, slopesDispersionInfoCtlSlim)

#separate out convergence vs divergence datasets
slopesDispersionRRconverge <- subset(slopesDispersionInfoRR, subset=(slope<0))
slopesDispersionRRconverge$converge_diverge <- 'converge'
slopesDispersionRRdiverge <- subset(slopesDispersionInfoRR, subset=(slope>0))
slopesDispersionRRdiverge$converge_diverge <- 'diverge'

#get response ratios
slopesDispersionRRconverge$RR <- with(slopesDispersionRRconverge, ifelse(ctl_slope>0, -1*(slope - ctl_slope) / ctl_slope,                                                                      (slope - ctl_slope) / ctl_slope))
slopesDispersionRRdiverge$RR <- with(slopesDispersionRRdiverge, ifelse(ctl_slope<0, -1*(slope - ctl_slope) / ctl_slope,
                                                           (slope - ctl_slope) / ctl_slope))

#merge experiment info with the response ratio dataframes
slopesDispersionInfoRRconverge <- merge(slopesDispersionRRconverge, info)
slopesDispersionInfoRRdiverge <- merge(slopesDispersionRRdiverge, info)

#number of factors model
summary(lm(RR ~ plot_mani.x, data=slopesDispersionInfoRRconverge))
ggplot(data=slopesDispersionInfoRRconverge, aes(x=plot_mani.x, y=RR)) +
  geom_point()

summary(lm(RR ~ plot_mani.x, data=slopesDispersionInfoRRdiverge))
ggplot(data=slopesDispersionInfoRRdiverge, aes(x=plot_mani.x, y=RR)) +
  geom_point()

#gamma diversity model
summary(lm(RR ~ species_num, data=slopesDispersionInfoRRconverge))
ggplot(data=slopesDispersionInfoRRconverge, aes(x=species_num, y=RR)) +
  geom_point()

summary(lm(RR ~ species_num, data=slopesDispersionInfoRRdiverge))
ggplot(data=slopesDispersionInfoRRdiverge, aes(x=species_num, y=RR)) +
  geom_point()

#MAP model
summary(lm(RR ~ MAP, data=slopesDispersionInfoRRconverge))
ggplot(data=slopesDispersionInfoRRconverge, aes(x=MAP, y=RR)) +
  geom_point()

summary(lm(RR ~ MAP, data=slopesDispersionInfoRRdiverge))
ggplot(data=slopesDispersionInfoRRdiverge, aes(x=MAP, y=RR)) +
  geom_point()

#ANPP model
summary(lm(RR ~ ANPP, data=slopesDispersionInfoRRconverge))
ggplot(data=slopesDispersionInfoRRconverge, aes(x=ANPP, y=RR)) +
  geom_point()

summary(lm(RR ~ ANPP, data=slopesDispersionInfoRRdiverge))
ggplot(data=slopesDispersionInfoRRdiverge, aes(x=ANPP, y=RR)) +
  geom_point()












