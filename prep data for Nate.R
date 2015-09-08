library(plyr)

setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

#read in experiment info
info <- read.csv('exp_info072015.csv')

#read in full change in mean and dispersion dataset
all <- read.csv('dispersion_and_means_press_experiments_with_exp_info_08062015.csv')

#get avergage dispersion
keep <- c('label', 'expt.year', 'dist', 'treatment', 'cal.year', 'trt.year', 'expt', 'community_type', 'project_name', 'site_code', 'data_type', 'MAP', 'ANPP', 'broad.ecosystem.type', 'ecosystem.type', 'species_num', 'nutrients', 'light', 'water', 'carbon', 'other_manipulation', 'num_manipulations')
disp <- subset(all, mean.disp=='disp')
disp <- disp[,colnames(disp) %in% keep]
dispMean <- ddply(disp, c('label', 'expt.year', 'treatment', 'cal.year', 'trt.year',  'expt', 'community_type', 'project_name', 'site_code', 'data_type', 'MAP', 'ANPP', 'broad.ecosystem.type', 'ecosystem.type', 'species_num', 'nutrients', 'light', 'water', 'carbon', 'other_manipulation', 'num_manipulations'), summarise,
                  dispersion=mean(dist))

#get mean data
mean <- subset(all, mean.disp=='mean')

#merge together mean and dispersion data
allClean <- merge(dispMean, mean, by=c('label', 'expt.year', 'treatment', 'cal.year', 'trt.year', 'expt', 'community_type', 'project_name', 'site_code', 'data_type', 'MAP', 'ANPP', 'broad.ecosystem.type', 'ecosystem.type', 'species_num', 'nutrients', 'light', 'water', 'carbon', 'other_manipulation', 'num_manipulations'), all=T)
names(allClean)[names(allClean)=="dist"] <- "mean_change"
names(allClean)[names(allClean)=="plot_mani.x"] <- "plot_mani"
drop <- c('X', 'plot_id', 'plot_mani.y', 'dist.1', 'dist.2', 'mean.disp', 'dg', 'site', 'dataset_site_code', 'true_num_manipulations', 'true_plot_mani')
allClean <- allClean[,!colnames(allClean) %in% drop]

#export to dropbox
write.csv(allClean, 'dispersion_means_press_bayesian_08172015.csv')





