library(plyr)
library(vegan)
library(lme4)
library(ggplot2)

setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank())

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
###################################################################################################
###################################################################################################
###################################################################################################


#read in experiment info
info <- read.csv('exp_info072015.csv')

#read in full spp cover dataset
all <- read.csv('all_relcov2_08062015.csv')

#create species matrix with one identifier per plot
drop <- c('X', 'unid', 'id', 'site_code', 'project_name', 'nutrients', 'light', 'carbon', 'water', 'other_manipulation', 'num_manipulations', 'clip', 'temp', 'precip', 'plot_mani', 'calendar_year', 'experiment_year', 'treatment_year', 'block', 'treatment', 'plot_id', 'plot_id1', 'data_type', 'species_num', 'n', 'p', 'k', 'herb_removal', 'burn', 'true_num_manipulations', 'c', 'plant_mani', 'true_plot_mani', 'community_type', 'lime', 'other_nut', 'cessation', 'dist', 'precip_vari', 'precip_vari_season', 'patchiness', 'other_manipulations', 'l', 'fungicide', 'soil_carbon', 'grazed', 'soil_depth', 'precip_season')
sppMatrix <- all[,!colnames(all) %in% drop]
rowInfo <- all[,colnames(all) %in% drop]

#calculate diversity, richness, evenness (Pielou's)
H <- diversity(sppMatrix)
S <- specnumber(sppMatrix)
J <- H/log(S)
diversity <- cbind(rowInfo, H, S, J)
diversity$J <- ifelse(diversity$J=='NaN', 1, diversity$J) #for plots with just one species...  should this be 0 or 1, or just dropped?

#calculate response ratios
diversity$expt <- with(diversity, paste(site_code, project_name, sep="_"))
meanDiversity <- ddply(diversity, c('expt', 'treatment_year', 'treatment', 'plot_mani'), summarise,
                       S=mean(S), 
                       J=mean(J))
trtMeanDiversity <- subset(meanDiversity, subset=(plot_mani>0))
ctlMeanDiversity <- subset(meanDiversity, subset=(plot_mani==0))
drop <- c('treatment', 'plot_mani')
ctlMeanDiversity <- ctlMeanDiversity[,!colnames(ctlMeanDiversity) %in% drop]
names(ctlMeanDiversity)[names(ctlMeanDiversity)=='S'] <- 'S_ctl'
names(ctlMeanDiversity)[names(ctlMeanDiversity)=='J'] <- 'J_ctl'
rrDiversity <- merge(trtMeanDiversity, ctlMeanDiversity, by=c('expt', 'treatment_year'))
rrDiversity$S_lnRR <- with(rrDiversity, log(S/S_ctl))
rrDiversity$J_lnRR <- with(rrDiversity, log(J/J_ctl))

#mixed effects model for richness and evenness 
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessModel <- lmer(S_lnRR ~ treatment_year*plot_mani + (1|expt), data=subset(rrDiversity, treatment_year>0)))
summary(evennessModel <- lmer(J_lnRR ~ treatment_year*plot_mani + (1|expt), data=subset(rrDiversity, treatment_year>0)))

#plot richness responses
colorManiDiverge <- c('#313695', '#4575b4', '#74add1', '#abd9e9', '#f46d43', '#d73027', '#a50026')
ggplot(data=subset(rrDiversity, treatment_year>0), aes(x=treatment_year, y=S_lnRR)) +
  geom_point(aes(y=S_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani)) +
  scale_x_continuous('Treatment Year')
ggplot(data=subset(rrDiversity, treatment_year>0), aes(x=factor(plot_mani), y=S_lnRR)) +
  geom_boxplot() +
  xlab('Number of Factors Manipulated')

#plot evenness responses
ggplot(data=subset(rrDiversity, treatment_year>0), aes(x=treatment_year, y=J_lnRR)) +
  geom_point(aes(y=J_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani)) +
  geom_smooth(aes(y=J_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani), method=lm, formula=y~log(x), se=T, size=3)
  scale_x_continuous('Treatment Year')

