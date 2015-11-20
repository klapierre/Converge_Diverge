setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(vegan)

data<-read.csv("all_relcov2_08062015.csv")
expinfo<-read.csv("exp_info072015.csv")
species<-read.csv("AllSpDataRaw_11192015.csv")

exp_info<-data%>%
 select(site_code, project_name, community_type, calendar_year, nutrients, light, carbon, water, other_manipulation, num_manipulations, clip, temp, precip, plot_mani, treatment, data_type, n, p, k, herb_removal, burn, true_num_manipulations, c, plant_mani, true_plot_mani, lime, other_nut, cessation, dist, precip_vari, precip_vari_season, patchiness, other_manipulations, l, fungicide, soil_carbon, grazed, soil_depth, precip_season, -X)%>%
  unique()

e001<-data%>%
  filter(project_name=="e001")%>%
  select(community_type, treatment)%>%
  unique

#write this file and manually add in MNR, RIO, CDR-e001

write.csv(exp_info, "experiment_information.csv")
read.csv("experiment_information2.csv")

experiment_length<-data%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year)%>%
  filter(treatment_year!=0)%>%
  unique%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(experiment_length=max(treatment_year))

species_num<-species%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, genus_species)%>%
  summarize(rich=length(abundance))%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(rich=length(rich))

merged<-merge(experiment_length, species_num, by=c("site_code","project_name","community_type"))
with(merged, cor.test(experiment_length, rich))# yes, experiment_legnth correlates with richness.
with(merged, plot(experiment_length, rich))


  
write.csv(exp_info, "Test.csv")

coltest<-exp_info%>%
  mutate(present=1)%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(presnt=length(present))

exp_info2<-merge(exp_info, expinfo, by=c("site_code","project_name","community_type"), all=T)
write.csv(exp_info2,'ExperimentInfo.csv')

coltest2<-exp_info2%>%
  mutate(present=1)%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(presnt2=length(present))

test<-merge(coltest, coltest2, by=c("site_code","project_name","community_type"))
#problems: in the source file alp_nme and lg_herbwood are listed twice
#the merge doubles ALP_NME for some reason, this take care of 8 of the extra reps.
#lg herbwood doubles for some reason in the merge, this is the extra 5 reps.
