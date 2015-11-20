library(plyr)
library(vegan)
library(lme4)
library(ggplot2)
library(RColorBrewer)
library(mgcv)

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

#read in change in mean and dispersion data
converge <- read.csv('dispersion_means_press_bayesian_08172015.csv')

#create species matrix with one identifier per plot
all$expt <- with(all, paste(site_code, project_name, community_type, sep="::"))
all$label <- with(all, paste(expt, treatment, sep="::"))
names(all)[names(all)=='treatment_year'] <- 'trt.year'
drop <- c('X', 'unid', 'id', 'site_code', 'project_name', 'nutrients', 'light', 'carbon', 'water', 'other_manipulation', 'num_manipulations', 'clip', 'temp', 'precip', 'plot_mani', 'calendar_year', 'experiment_year', 'block', 'treatment', 'plot_id', 'plot_id1', 'data_type', 'species_num', 'n', 'p', 'k', 'herb_removal', 'burn', 'true_num_manipulations', 'c', 'plant_mani', 'true_plot_mani', 'community_type', 'lime', 'other_nut', 'cessation', 'dist', 'precip_vari', 'precip_vari_season', 'patchiness', 'other_manipulations', 'l', 'fungicide', 'soil_carbon', 'grazed', 'soil_depth', 'precip_season', 'expt', 'label', 'trt.year')
sppMatrix <- all[,!colnames(all) %in% drop]
rowInfo <- all[,colnames(all) %in% drop]

#calculate diversity, richness, evenness (Pielou's)
H <- diversity(sppMatrix)
S <- specnumber(sppMatrix)
J <- H/log(S)
diversity <- cbind(rowInfo, H, S, J)
diversity$J <- ifelse(diversity$J=='NaN', 0, diversity$J) #for plots with just one species...  should this be 0 or 1, or just dropped?

#merge with change in mean and dispersion data
meanDiversity <- ddply(diversity, c('expt', 'label', 'trt.year', 'treatment', 'plot_mani'), summarise,
                       S_mean=mean(S), 
                       J_mean=mean(J),
                       S_sd=sd(S),
                       J_sd=sd(J),
                       S_N=length(S),
                       J_N=length(J))
keep <- c('expt', 'trt.year', 'label', 'treatment', 'J_mean', 'S_mean', 'S_sd', 'J_sd', 'S_N', 'J_N')
meanDiversityMerge <- meanDiversity[,colnames(meanDiversity) %in% keep]
meanDiversityMerge2 <- merge(meanDiversityMerge, converge, by=c('expt', 'label', 'treatment', 'trt.year'))
write.csv(meanDiversityMerge2, 'dispersion_means_press_bayesian_09092015.csv')

#calculate response ratios
trtMeanDiversity <- subset(meanDiversity, subset=(plot_mani>0))
ctlMeanDiversity <- subset(meanDiversity, subset=(plot_mani==0))
drop <- c('treatment', 'plot_mani')
ctlMeanDiversity <- ctlMeanDiversity[,!colnames(ctlMeanDiversity) %in% drop]
names(ctlMeanDiversity)[names(ctlMeanDiversity)=='S_mean'] <- 'S_ctl'
names(ctlMeanDiversity)[names(ctlMeanDiversity)=='J_mean'] <- 'J_ctl'
names(ctlMeanDiversity)[names(ctlMeanDiversity)=='S_sd'] <- 'S_sd_ctl'
names(ctlMeanDiversity)[names(ctlMeanDiversity)=='J_sd'] <- 'J_sd_ctl'
rrDiversity <- merge(trtMeanDiversity, ctlMeanDiversity, by=c('expt', 'trt.year'))
rrDiversity$S_lnRR <- with(rrDiversity, log(S_mean/S_ctl))
rrDiversity$J_lnRR <- with(rrDiversity, log(J_mean/J_ctl))
rrDiversity$S_percent <- with(rrDiversity, 100*(S_mean-S_ctl)/S_ctl)
rrDiversity$J_percent <- with(rrDiversity, 100*(J_mean-J_ctl)/J_ctl)
rrDiversity$S_sd_percent <- with(rrDiversity, 100*(S_sd-S_sd_ctl)/S_sd_ctl)
rrDiversity$J_sd_percent <- with(rrDiversity, 100*(J_sd-J_sd_ctl)/J_sd_ctl)

###mixed effects model for richness and evenness - ln response ratio
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessModel <- lmer(S_lnRR ~ trt.year*plot_mani + (1|expt), data=subset(rrDiversity, trt.year>0)))
summary(evennessModel <- lmer(J_lnRR ~ trt.year*plot_mani + (1|expt), data=subset(rrDiversity, trt.year>0)))

#plot richness responses
colorManiDiverge <- c('#313695', '#4575b4', '#74add1', '#abd9e9', '#f46d43', '#d73027', '#a50026')
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=trt.year, y=S_lnRR)) +
  geom_point(aes(y=S_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani)) +
  #   geom_smooth(aes(y=S_lnRR, colour=plot_mani, group=plot_mani), method=lm, formula=y~log(x), se=T, size=3)
  xlab('Treatment Year') +
  ylab('ln Response Ratio (Richness)')
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=factor(plot_mani), y=S_lnRR)) +
  geom_boxplot() +
  xlab('Number of Factors Manipulated') +
  ylab('ln Response Ratio (Richness)')
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=plot_mani, y=S_lnRR)) +
  geom_point() +
  geom_smooth(method=lm, formula=y~log(x)) +
  xlab('Number of Factors Manipulated') +
  ylab('ln Response Ratio (Richness)')

#plot evenness responses
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=trt.year, y=J_lnRR)) +
  geom_point(aes(y=J_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani)) +
  geom_smooth(aes(y=J_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani), method=lm, formula=y~log(x), se=T, size=3)
scale_x_continuous('Treatment Year')


###mixed effects model for richness and evenness - percent change
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessPercentModel <- lmer(S_percent ~ trt.year*plot_mani + (1|expt), data=subset(rrDiversity, trt.year>0)))
summary(evennessPercentModel <- lmer(J_percent ~ trt.year*plot_mani + (1|expt), data=subset(rrDiversity, trt.year>0)))

#plot richness responses
rrDiversity$num_mani <- factor(as.character(rrDiversity$plot_mani),levels=as.character(c(1:7)))
colorManiDiverge <- c('#313695', '#4575b4', '#74add1', '#abd9e9', '#f46d43', '#d73027', '#a50026')
cg <- c(brewer.pal(7,"RdBu"))[7:1]
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=trt.year, y=S_percent)) +
  stat_smooth(aes(y=S_percent, group=interaction(expt,num_mani), colour=num_mani), method=lm, formula=y~log(x), se=F, alpha=0.7) +
  stat_smooth(aes(y=S_percent, group=num_mani, colour=num_mani, fill=num_mani), method='gam', formula=y~s(x, k=3), se=T, alpha=0.7) +
  #geom_ribbon(aes(y=S_percent, ymin=S_percent-2*S_sd_percent, ymax=S_percent+2*S_sd_percent, fill=num_mani))
  scale_colour_manual(values=cg, name="Number of manipulations") +
  scale_fill_manual(values=cg, name="Number of manipulations") +
  xlab('Treatment Year') +
  ylab('Change in Richness (%)')
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=factor(plot_mani), y=S_percent)) +
  geom_boxplot() +
  xlab('Number of Factors Manipulated') +
  ylab('Change in Richness (%)')
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=plot_mani, y=S_percent)) +
  geom_point() +
  geom_smooth(method=lm, formula=y~log(x)) +
  xlab('Number of Factors Manipulated') +
  ylab('Change in Richness (%)')

#plot evenness responses
ggplot(data=subset(rrDiversity, trt.year>0), aes(x=trt.year, y=J_lnRR)) +
  geom_point(aes(y=J_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani)) +
  geom_smooth(aes(y=J_lnRR, colour=plot_mani, fill=plot_mani, group=plot_mani), method=lm, formula=y~log(x), se=T, size=3)
scale_x_continuous('Treatment Year')

