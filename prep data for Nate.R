library(plyr)
library(vegan)
library(lme4)
library(ggplot2)
library(RColorBrewer)
library(mgcv)

setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

#read in experiment info
info <- read.csv('exp_info072015.csv')

#read in full change in mean and dispersion dataset
multivariate <- read.csv('dispersion_and_means_press_experiments_with_exp_info_08062015.csv')

#read in full spp cover dataset
all <- read.csv('all_relcov2_08062015.csv')

#get avergage dispersion
keep <- c('label', 'expt.year', 'dist', 'treatment', 'cal.year', 'trt.year', 'expt', 'community_type', 'project_name', 'site_code', 'data_type', 'MAP', 'ANPP', 'broad.ecosystem.type', 'ecosystem.type', 'species_num', 'nutrients', 'light', 'water', 'carbon', 'other_manipulation', 'num_manipulations')
disp <- subset(multivariate, mean.disp=='disp')
disp <- disp[,colnames(disp) %in% keep]
dispMean <- ddply(disp, c('label', 'expt.year', 'treatment', 'cal.year', 'trt.year',  'expt', 'community_type', 'project_name', 'site_code', 'data_type', 'MAP', 'ANPP', 'broad.ecosystem.type', 'ecosystem.type', 'species_num', 'nutrients', 'light', 'water', 'carbon', 'other_manipulation', 'num_manipulations'), summarise,
                  dispersion=mean(dist))

#get mean data
mean <- subset(multivariate, mean.disp=='mean')

#merge together mean and dispersion data
convergeClean <- merge(dispMean, mean, by=c('label', 'expt.year', 'treatment', 'cal.year', 'trt.year', 'expt', 'community_type', 'project_name', 'site_code', 'data_type', 'MAP', 'ANPP', 'broad.ecosystem.type', 'ecosystem.type', 'species_num', 'nutrients', 'light', 'water', 'carbon', 'other_manipulation', 'num_manipulations'), all=T)
names(convergeClean)[names(convergeClean)=="dist"] <- "mean_change"
names(convergeClean)[names(convergeClean)=="plot_mani.x"] <- "plot_mani"
drop <- c('X', 'plot_id', 'plot_mani.y', 'dist.1', 'dist.2', 'mean.disp', 'dg', 'site', 'dataset_site_code', 'true_num_manipulations', 'true_plot_mani')
convergeClean <- convergeClean[,!colnames(convergeClean) %in% drop]

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
invD <- diversity(sppMatrix, 'inv')
evenness <- invD/S
diversity <- cbind(rowInfo, H, S, evenness)
diversity$J <- ifelse(diversity$evenness=='NaN', 0, diversity$J) #for plots with just one species...  should this be 0 or 1, or just dropped?

#merge with change in mean and dispersion data
meanDiversity <- ddply(diversity, c('expt', 'label', 'trt.year', 'treatment', 'plot_mani'), summarise,
                       S_mean=mean(S), 
                       evenness_mean=mean(evenness),
                       H_mean=mean(H),
                       S_sd=sd(S),
                       evenness_sd=sd(evenness),
                       H_sd=sd(H),
                       S_N=length(S),
                       evenness_N=length(evenness),
                       H_N=length(H))
keep <- c('expt', 'trt.year', 'label', 'treatment', 'H_mean', 'evenness_mean', 'S_mean', 'H_sd', 'S_sd', 'evenness_sd', 'S_N', 'evenness_N', 'H_N')
meanDiversityMerge <- meanDiversity[,colnames(meanDiversity) %in% keep]
meanDiversityMerge2 <- merge(meanDiversityMerge, convergeClean, by=c('expt', 'label', 'treatment', 'trt.year'))

#subset out treatments where no resource was manipulated
allCleanResource <- subset(meanDiversityMerge2, subset=(label!='Alps::NME::0::CoG' & label!='ASGA::clonal::0::non-clonal_CO' & label!='BAY::LIND::0::ref_rich2' & label!='BAY::LIND::0::ref_rich4' & label!='BAY::LIND::0::ref_rich8' & label!='BAY::LIND::0::ref_rich16' & label!='BUX::PQ::0::warm' & label!='CEH::MEGARICH::0::AcEt' & label!='CLE::Imagine::0::T' & label!='CUL::Culardoch::0::clip' & label!='Finse::WarmNut::0::warming' & label!='JSP::GCE::0::H' & label!='KAEFS::WAPAClip::0::U WC' & label!='KAEFS::WAPAClip::0::C CC' & label!='KAEFS::WAPAClip::0::C WC' & label!='KLUB::BFFert::0::N0F1' & label!='KLUG::KGFert::0::N0B1' & label!='KNZ::BGP::0::u_m_c' & label!='KNZ::BGP::0::b_u_c' & label!='KNZ::BGP::0::b_m_c' & label!='KNZ::RaMPs::0::AH' & label!='KNZ::RHPs::0::stone' & label!='LATNJA::CLIP::Heath::T' & label!='maerc::fireplots::0::suuu' & label!='maerc::fireplots::0::wuuu' & label!='maerc::fireplots::0::uuug' & label!='maerc::fireplots::0::suug' & label!='maerc::fireplots::0::wuug' & label!='Manitoba::CCD::0::CHA' & label!='Manitoba::CCD::0::CLA' & label!='Manitoba::CCD::0::WHA' & label!='Manitoba::CCD::0::WLA' & label!='Manitoba::CCD::0::WNA' & label!='NGBER::gb::0::CURRENT' & label!='NIN::HerbDiv::0::2NF' & label!='NIN::HerbDiv::0::3NF' & label!='NIN::HerbDiv::0::4NF' & label!='NIN::HerbDiv::0::5NF' & label!='NWT::snow::0::XXW' & label!='Saskatchewan::CCD::0::CHA' & label!='Saskatchewan::CCD::0::CLA' & label!='Saskatchewan::CCD::0::WHA' & label!='Saskatchewan::CCD::0::WLA' & label!='Saskatchewan::CCD::0::WNA' & label!='SEV::WENNDEx::0::T' & label!='SKY::UK::0::H' & label!='Alberta::CCD::0::CHA' & label!='Alberta::CCD::0::CLA' & label!='Alberta::CCD::0::WHA' & label!='Alberta::CCD::0::WLA' & label!='Alberta::CCD::0::WNA' & trt.year!=0))

#export to dropbox
write.csv(allCleanResource, 'dispersion_means_press_bayesian_11172015.csv')





