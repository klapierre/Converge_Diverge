setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(vegan)

data<-read.csv("all_relcov2_08062015.csv")#old species data but has all the columns
siteinfo<-read.csv("SiteInfo_11202015.csv")%>%
  select(-X, -species_num)
species<-read.csv("SpeciesRawAbundance_11202015.csv")

exp_info<-data%>%
 select(site_code, project_name, community_type, calendar_year, treatment_year, nutrients, light, carbon, water, other_manipulation, num_manipulations, clip, temp, precip, plot_mani, treatment, data_type, n, p, k, herb_removal, burn, true_num_manipulations, c, plant_mani, true_plot_mani, lime, other_nut, cessation, dist, precip_vari, precip_vari_season, patchiness, other_manipulations, l, fungicide, soil_carbon, grazed, soil_depth, precip_season, -X)%>%
  unique()

mnr<-read.csv("~/Dropbox/converge_diverge/datasets/ORIGINAL_DATASETS/MNR_watfer/Exp_Info.csv")%>%
  mutate(clip=0, community_type=0, temp=0, herb_removal=0, burn=0, true_num_manipulations=0, c=0, plant_mani=0, true_plot_mani=0, lime=0, cessation=0, dist=0, precip_vari=0, precip_vari_season=0, patchiness=0, other_manipulations=0, l=0, fungicide=0, soil_carbon=0, grazed=0, soil_depth=0, precip_season=0)%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, nutrients, light, carbon, water, other_manipulation, num_manipulations, clip, temp, precip, plot_mani, treatment, data_type, n, p, k, herb_removal, burn, true_num_manipulations, c, plant_mani, true_plot_mani, lime, other_nut, cessation, dist, precip_vari, precip_vari_season, patchiness, other_manipulations, l, fungicide, soil_carbon, grazed, soil_depth, precip_season, -X)%>%
  unique()

rio<-read.delim("~/Dropbox/converge_diverge/datasets/FINAL_SEPT2014/clean datasets - please do not touch/sp text files/RIO_interaction.txt")%>%
  mutate(nutrients=1, light=0, carbon=0, water=1, other_manipulation=0, num_manipulations=2, clip=0, temp=0, herb_removal=0, burn=0, true_num_manipulations=0, c=0, plant_mani=0, true_plot_mani=0, lime=0, cessation=0, dist=0, precip_vari_season=0, p=0, k=0, other_nut=0, patchiness=0, other_manipulations=0, l=0, fungicide=0, soil_carbon=0, grazed=0, soil_depth=0, precip_season=0)%>%
  select(site_code, project_name, community_type, calendar_year, nutrients, treatment_year, light, carbon, water, other_manipulation, num_manipulations, clip, temp, precip, plot_mani, treatment, data_type, n, p, k, herb_removal, burn, true_num_manipulations, c, plant_mani, true_plot_mani, lime, other_nut, cessation, dist, precip_vari, precip_vari_season, patchiness, other_manipulations, l, fungicide, soil_carbon, grazed, soil_depth, precip_season)%>%
  unique()

exp_info2<-rbind(exp_info, rio, mnr)
write.csv(exp_info2, "ExperimentInformation_11202015.csv")
  
experiment_length<-exp_info2%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(experiment_length=max(treatment_year))

siteinfo2<-merge(siteinfo, experiment_length, by=c("site_code","project_name","community_type"), all=T)

species_num<-species%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, genus_species)%>%
  summarize(rich=length(abundance))%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(rich=length(rich))

site_info3<-merge(siteinfo2, species_num, by=c("site_code","project_name","community_type"), all=T)

##need to rarefy species. not sure how to do this with cover data.
