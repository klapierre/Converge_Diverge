library(ggplot2)
library(plyr)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

setwd('C:\\Users\\Kim\\Dropbox\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

#read in Forest's ASReml output
dispersion <- read.csv('WithinTrtDivergenceConvergenceRates.csv')

#read in Meghan's data info
info <- read.csv('exp_info072015.csv')

#read in full change in mean and dispersion dataset
full <- read.csv('dispersion_and_means_press_experiments_with_exp_info_03232015.csv')

#create column of number of factors manipulated (rough estimate based on each resource column)
dispersion$factor_num <- with(dispersion, nut2+pp+car+heat+other)

#get treatment and control effects as columns
dispersionTrt <- subset(dispersion, subset=(binary_trt!='0.0.0.0.0'))
dispersionCtl <- subset(dispersion, subset=(binary_trt=='0.0.0.0.0'))
dispersionCtlSlim <- dispersionCtl[,-c(2, 4:9)]
names(dispersionCtlSlim)[names(dispersionCtlSlim)=='eff'] <- 'ctl_eff'
dispersionRR <- merge(dispersionTrt, dispersionCtlSlim)

#separate out convergence vs divergence datasets
dispersionRRconverge <- subset(dispersionRR, subset=(eff<0))
dispersionRRconverge$converge_diverge <- 'converge'
dispersionRRdiverge <- subset(dispersionRR, subset=(eff>0))
dispersionRRdiverge$converge_diverge <- 'diverge'

#get response ratios
dispersionRRconverge$RR <- with(dispersionRRconverge, ifelse(ctl_eff>0, -1*(eff - ctl_eff) / ctl_eff,                                                                      (eff - ctl_eff) / ctl_eff))
dispersionRRdiverge$RR <- with(dispersionRRdiverge, ifelse(ctl_eff<0, -1*(eff - ctl_eff) / ctl_eff,
                                                           (eff - ctl_eff) / ctl_eff))

#recombinde converge and diverge response ratio dataframes
dispersionRRcombo <- rbind(dispersionRRconverge, dispersionRRdiverge)

#merge experiment info with the response ratio dataframe
dispersionRRinfo <- merge(dispersionRRcombo, info)

# #regress RR with gamma diversity
# gammaModelConverge <- lm(RR ~ species_num, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(gammaModelConverge)
# 
# gammaModelDiverge <- lm(RR ~ species_num, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(gammaModelDiverge)
# 
# ggplot(data=dispersionRRinfo, aes(x=species_num, y=RR, colour=converge_diverge)) +
#   geom_point(size=5)
# 
# 
# #regress RR with MAP
# MAPModelConverge <- lm(RR ~ MAP, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(MAPModelConverge)
# 
# MAPModelDiverge <- lm(RR ~ MAP, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(MAPModelDiverge)
# 
# ggplot(data=dispersionRRinfo, aes(x=MAP, y=RR, colour=converge_diverge)) +
#   geom_point(size=5)
# 
# 
# #regress RR with ANPP
# ANPPModelConverge <- lm(RR ~ ANPP, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(ANPPModelConverge)
# 
# ANPPModelDiverge <- lm(RR ~ ANPP, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(ANPPModelDiverge)
# 
# ggplot(data=dispersionRRinfo, aes(x=ANPP, y=RR, colour=converge_diverge)) +
#   geom_point(size=5)
# 
# 
# #regress RR with number of factors manipulated
# factorsModelConverge <- lm(RR ~ factor_num, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(factorsModelConverge)
# 
# factorsModelDiverge <- lm(RR ~ factor_num, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(factorsModelDiverge)
# 
# ggplot(data=dispersionRRinfo, aes(x=factor_num, y=RR, colour=converge_diverge)) +
#   geom_point(size=5)
# 
# ###multiple regression with number of factors included
# #gamma diversity and number of factors
# gammaFactorModelConverge <- lm(RR ~ factor_num*species_num, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(gammaFactorModelConverge)
# 
# gammaFactorModelDiverge <- lm(RR ~ factor_num*species_num, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(gammaFactorModelDiverge)
# 
# #ANPP and number of factors
# ANPPFactorModelConverge <- lm(RR ~ factor_num*ANPP, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(ANPPFactorModelConverge)
# 
# ANPPFactorModelDiverge <- lm(RR ~ factor_num*ANPP, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(ANPPFactorModelDiverge)
# 
# #MAP and number of factors
# MAPFactorModelConverge <- lm(RR ~ factor_num*MAP, data=subset(dispersionRRinfo, converge_diverge=='converge'))
# summary(MAPFactorModelConverge)
# 
# ggplot(data=subset(dispersionRRinfo, converge_diverge=='converge'), aes(x=MAP, y=RR, colour=factor_num)) +
#   geom_point(size=5)
# 
# MAPFactorModelDiverge <- lm(RR ~ factor_num*MAP, data=subset(dispersionRRinfo, converge_diverge=='diverge'))
# summary(MAPFactorModelDiverge)

#make dataset indicating what control treatments are, plus all other relevant information
controls <- full[,-c(1:5, 7:8, 10:74)]
controlsAggregate <- aggregate(plot_mani.x ~ label, mean, data=full)

#using full dataset, get slopes of dispersion through time for each treatment independantly
slopesDispersion <- ddply(subset(full, mean.disp=='disp'), 'label', function(x) {
  model <- lm(dist~trt.year, data=x)
  coef(model)
names(slopesDispersion)[names(slopesDispersion)=="trt.year"] <- "slope"

#merge slopes data with controls and extra info
slopesDispersionInfo <- merge(slopesDispersion, controlsAggregate)












