#Kim's
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\clean datasets - please do not touch\\sp text files')

#Meghan's
setwd("~/Dropbox/converge_diverge/datasets/FINAL_SEPT2014/clean datasets - please do not touch/sp text files")

library(tidyr)
library(dplyr)

# notes:
# nutrients, light, carbon, water, and other are binary variables for the entire experiment (whether one of these factors was manipulated)
# n-other_trt are variables describing the specific treatment amounts or categories
#### if all plots were burned (at any frequency (<20yrs) over the course of the experiment or at the site where the experiment takes place), fenced to remove herbivores (at any point in time where the experiment takes place), or mowed/clipped, then the variable gets a 1 for all plots. THESE are also treatments for some experiments.
#### successional and plant_mani are binary variables
####   if successional=1, then all plots were disturbed in some way prior to the experiment start
#   if plant_mani=1, then some or all species are planted into the plots at the start of the experiment
#### cessation is a binary variable for each treatment and year combination, and is 1 when the treatment has been stopped - DROPPING CESSATION
# plot_mani is the total number of factors manipulated compared to the control (i.e., if all plots are burned plot_mani does not increase, but if treatment plots are burned when controls are not, then plot_mani increases by 1)
# resource_mani is a binary variable to subset out the treatments that do not directly manipulate a resource (e.g., only increased temperature)
#   resource_mani=1 for all control plots and any treatment plot that directly manipulates a resource; resource_mani=0 for plots where no resource was directly manipulated (except for the controls)
# pulse is a binary variable for treatments that were a one time application that did not also get applied to the control plots (e.g., burning in CUL, which was a treatment just in the first year of the experiment, but was not applied to the controls)
#   a factor was not considered a pulse if all plots experienced it (e.g., solarization in ASGA clonal, wildfire in SEV WENNDex)
# public means the dataset is publically available and we don't need to ask permission to use it

watering<-read.delim("ANG_watering.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='W', 20, ifelse(treatment=='S', 20, 0)), 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=ifelse(treatment=='W', 'winter water addition', ifelse(treatment=='S', 'spring water addition', 0)), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='W', 1, ifelse(treatment=='S', 1, 0)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

mat2<-read.delim("ARC_mat2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NP', 10, 0), 
         p=ifelse(treatment=='NP', 5, 0), 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='NP', 2, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

mnt<-read.delim("ARC_mnt.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NP', 10, 0), 
         p=ifelse(treatment=='NP', 5, 0), 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='NP', 2, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

clonal<-read.delim("ASGA_Clonal.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='mixed_CO', 0, ifelse(treatment=='non-clonal_CO', 0, 20.1)), 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0,
         trt_details=ifelse(treatment=='non-clonal_CO', 'non-clonal species', ifelse(treatment=='mixed_LP', 'large nutrient patches', ifelse(treatment=='non-clonal_LP', 'non-clonal species, large nutrient patches', ifelse(treatment=='mixed_SP', 'small nutrient patches', ifelse(treatment=='non-clonal_SP', 'non-clonal species, small nutrient patches', ifelse(treatment=='non-clonal_UN', 'non-clonal species', 0)))))),
         successional=1, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='non-clonal_UN'|treatment=='non-clonal_LP'|treatment=='non-clonal_SP', 2, ifelse(treatment=='non-clonal_CO'|treatment=='mixed_LP'|treatment=='mixed_SP'|treatment=='mixed_UN', 1, 0)))%>%
  mutate(resource_mani=ifelse(treatment=='non-clonal_CO', 0, 1))%>%
  mutate(public=0)%>%
  unique()

exp1<-read.delim("ASGA_Exp1.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='2_0_CO'|treatment=='1_0_CO'|treatment=='2_1_CO'|treatment=='1_1_CO', 0, 20.1),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         other_trt=0,
         trt_details=ifelse(treatment=='2_0_PA'|treatment=='1_0_PA'|treatment=='2_1_PA'|treatment=='1_1_PA', 'nutrient patches', 0),
         successional=ifelse(treatment=='2_0_CO'|treatment=='2_0_PA'|treatment=='2_0_UN'|treatment=='2_1_CO'|treatment=='2_1_PA'|treatment=='2_1_UN', 1, 0),
         plant_mani=ifelse(treatment=='2_1_CO'|treatment=='2_1_PA'|treatment=='2_1_UN'|treatment=='1_1_CO'|treatment=='1_1_PA'|treatment=='1_1_UN', 1, 0),
         pulse=ifelse(treatment=='2_0_PA'|treatment=='2_0_UN'|treatment=='2_0_CO'|treatment=='2_1_PA'|treatment=='2_1_UN'|treatment=='2_1_CO', 1, 0))%>%
  mutate(plot_mani=ifelse(treatment=='1_0_CO', 0, ifelse(treatment=='1_0_PA'|treatment=='1_0_UN'|treatment=='1_1_CO'|treatment=='2_0_CO', 1, ifelse(treatment=='1_1_PA'|treatment=='1_1_UN'|treatment=='2_0_PA'|treatment=='2_0_UN'|treatment=='2_1_CO', 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='2_0_CO'|treatment=='1_1_CO'|treatment=='2_1_CO', 0, 1))%>%
  mutate(public=0)%>%
  unique()

nitphos<-read.csv("AZI_NitPhos.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0, 
         n=ifelse(treatment=='N1P0', 5, ifelse(treatment=='N2P0', 10, ifelse(treatment=="N3P0", 15, ifelse(treatment=="N0P0", 0, 10)))), 
         p=ifelse(treatment=='N2P1', 2, ifelse(treatment=="N2P2", 4, ifelse(treatment=="N2P3", 8, 0))), 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0P0', 0, ifelse(treatment=="N1P0", 1, ifelse(treatment=="N2P0", 1, ifelse(treatment=="N3P0", 1, 2)))))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

#16 spp plots are controls
lind<-read.delim("BAY_LIND.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='rain_rich1'|treatment=='rain_rich2'|treatment=='rain_rich4'|treatment=='rain_rich8'|treatment=='rain_rich16', 8, 0),
         temp=0, 
         mow_clip=1, 
         burn=0, 
         herb_removal=0, 
         pulse=0,
         other_trt=0,
         trt_details=ifelse(treatment=='ref_rich1'|treatment=='rain_rich1', '1 sp', ifelse(treatment=='ref_rich2'|treatment=='rain_rich2', '2 sp', ifelse(treatment=='ref_rich4'|treatment=='rain_rich4', '4 sp', ifelse(treatment=='ref_rich8'|treatment=='rain_rich8', '8 sp', '16 sp')))),
         successional=1, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='ref_rich16', 0, ifelse(treatment=='ref_rich1'|treatment=='ref_rich2'|treatment=='ref_rich4'|treatment=='ref_rich8', 1, 2)))%>%
  mutate(resource_mani=ifelse(treatment=='ref_rich1'|treatment=='ref_rich2'|treatment=='ref_rich4'|treatment=='ref_rich8', 0, 1))%>%
  mutate(public=0)%>%
  unique()

events<-read.delim("Bt_EVENT2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=1, 
         burn=0, 
         herb_removal=0,
         other_trt=ifelse(treatment=='CM-N1', 'reduced precip variability', ifelse(treatment=='D1-N1', 'early drought', ifelse(treatment=='D2-N1', 'late drought', 0))), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CA-N1', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

pq<-read.delim("BUX_PQ.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=='wet'|treatment=='warm wet', 20, ifelse(treatment=='dry'|treatment=='warm dry', -20, 0)),
         temp=ifelse(treatment=='warm'|treatment=='warm dry'|treatment=='warm wet', 3, 0),
         mow_clip=1, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='warm'|treatment=='dry'|treatment=='wet', 1, ifelse(treatment=='control', 0, 2)))%>%
  mutate(resource_mani=ifelse(treatment=='warm', 0, 1))%>%
  mutate(public=1)%>%
  unique()

pennings<-read.delim("CAR_Pennings.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
        n=ifelse(treatment=='NPK'&calendar_year>1999, 164, ifelse(treatment=='NPK'&calendar_year==1999, 84, 0)),
        p=ifelse(treatment=='NPK'&calendar_year>1999, 82, ifelse(treatment=='NPK'&calendar_year==1999, 42, 0)),
        k=ifelse(treatment=='NPK'&calendar_year>1999, 41, ifelse(treatment=='NPK'&calendar_year==1999, 21, 0)),
        CO2=0, 
        precip=0,
        temp=0,
        mow_clip=0, 
        burn=0, 
        herb_removal=0,
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0, 
        pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='NPK', 3, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

rmapc<-read.delim("CAU_RMAPC.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N', 9, ifelse(treatment=='NP', 9, 0)),
         p=ifelse(treatment=='P', 2.6, ifelse(treatment=='NP', 2.6, 0)),
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='H2O', 20, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=ifelse(treatment=='Ca', "lime added", 0), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Cont', 0, ifelse(treatment=='NP', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

biocon<-read.delim("CDR_biocon.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=1, water=0, other_manipulation=0,
         n=ifelse(treatment=='N_X', 4, ifelse(treatment=='N_C', 4, 0)),
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='X_C', 160, ifelse(treatment=='N_C', 160, 0)),
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0,
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='X_X', 0, ifelse(treatment=='N_C', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

e001<-read.csv("CDR_e001.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='2_y_n', 1.02, ifelse(treatment=='3_y_n', 2.04, ifelse(treatment=='4_y_n', 3.40, ifelse(treatment=='5_y_n', 5.44, ifelse(treatment=='6_y_n', 9.52, ifelse(treatment=='7_y_n', 17, ifelse(treatment=='8_y_n', 27.2, 0))))))),
         p=ifelse(treatment=='9_y_n', 0, 4.6), 
         k=ifelse(treatment=='9_y_n', 0, 6.1),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=1,
         other_trt=ifelse(treatment=='9_y_n', 0, "mirconutrients and lime added"), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='1_y_n', 4, ifelse(treatment=='9_y_n', 0, 5)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

e002<-read.delim("CDR_e002.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='2_f_u_n', 1.02, ifelse(treatment=='3_f_u_n', 2.04, ifelse(treatment=='4_f_u_n', 3.4, ifelse(treatment=='5_f_u_n', 5.44, ifelse(treatment=='6_f_u_n', 9.52, ifelse(treatment=='7_f_u_n', 17, ifelse(treatment=='8_f_u_n', 27.2, 0))))))),
         p=ifelse(treatment=='9_f_u_n', 0, 4.6),
         k=ifelse(treatment=='9_f_u_n', 0, 6.1),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=1,
         other_trt=ifelse(treatment=='9_f_u_n', 0, "micronutrients and lime added"),
         trt_details=0,
         successional=1, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='1_f_u_n', 4, ifelse(treatment=='9_f_u_n', 0, 5)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  filter(calendar_year<1992)%>%##drops everything once cessation starts
  unique()

megarich<-read.delim("CEH_Megarich.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
        nutrients=1, light=0, carbon=1, water=0, other_manipulation=1,
        n=10, 
        p=2, 
        k=2,
        CO2=ifelse(treatment=='EcAt', 280, ifelse(treatment=='EcEt', 280, 0)),
        precip=0,
        temp=ifelse(treatment=='AcEt', 2.9, ifelse(treatment=='EcEt', 2.9, 0)),
        mow_clip=1, 
        burn=0, 
        herb_removal=0,
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0, 
        pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='AcAt', 0, ifelse(treatment=='EcEt', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='AcEt', 0, 1))%>%
  mutate(public=0)%>%
  unique()

imagine<-read.delim("CLE_imagine.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=1, other_manipulation=1,
         n=0,
         p=0, 
         k=0,
         CO2=ifelse(treatment=='TDCO2', 200, 0),
         precip=ifelse(treatment=='TD', -20, ifelse(treatment=='TDCO2', -20, 0)),
         temp=ifelse(treatment=='C', 0, 3.5),
         mow_clip=1, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='T', 1, ifelse(treatment=='TD', 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='T', 0, 1))%>%
  mutate(public=0)%>%
  unique()

culardoch<-read.delim("CUL_culardoch.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N10'|treatment=='N10burn'|treatment=='N10clip'|treatment=='N10burnclip', 1, ifelse(treatment=='N20'|treatment=='N20burn'|treatment=='N20clip'|treatment=='N20burnclip', 2, ifelse(treatment=='N50'|treatment=='N50burn'|treatment=='N50clip'|treatment=='N50burnclip', 5, 0))),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment=='clip'|treatment=='burnclip'|treatment=='N10clip'|treatment=='N20clip'|treatment=='N50clip'|treatment=='N10burnclip'|treatment=='N20burnclip'|treatment=='N50burnclip', 1, 0),
         burn=ifelse(treatment=='N10burn'|treatment=='N20burn'|treatment=='N50burn'|treatment=='burn'|treatment=='burnclip'|treatment=='N10burnclip'|treatment=='N20burnclip'|treatment=='N50burnclip', 1, 0),
        herb_removal=0,
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0,
        pulse=ifelse(treatment=='burn'|treatment=='burnclip'|treatment=='N10burn'|treatment=='N10burnclip'|treatment=='N20burn'|treatment=='N20burnclip'|treatment=='N50burn'|treatment=='N50burnclip', 1, 0))%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment=='N10'|treatment=='N20'|treatment=='N50'|treatment=='burn'|treatment=='clip', 1, ifelse(treatment=='N10burnclip'|treatment=='N20burnclip'|treatment=='N50burnclip', 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment=='burn'|treatment=='clip'|treatment=='burnclip', 0, 1))%>%
  mutate(public=0)%>%
  unique()

gap2<-read.delim("DCGS_gap2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=1, carbon=0, water=0, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0,
         burn=0, 
         herb_removal=0, 
         other_trt=ifelse(treatment=='_018', '18 ft opening', ifelse(treatment=='_033', '33 ft opening', ifelse(treatment=='_066', '66 ft opening', ifelse(treatment=='_100', '100 ft opening', ifelse(treatment=='_150', '150 ft opening', 0))))),
         trt_details=0,
         successional=0, 
         plant_mani=0,
         pulse=1)%>%
  mutate(plot_mani=ifelse(treatment=='_000', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

nsfc<-read.delim("DL_NSFC.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N', 10, ifelse(treatment=='WN', 10, 0)),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='W', 49.8, ifelse(treatment=='WN', 49.8, 0)),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='WN', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

warmnut<-read.delim("Finse_WarmNut.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='nutrient addition'|treatment=='warming + nutrient addition', 10, 0),
         p=ifelse(treatment=='nutrient addition'|treatment=='warming + nutrient addition', 2, 0),
         k=ifelse(treatment=='nutrient addition'|treatment=='warming + nutrient addition', 8, 0),
         CO2=0, 
         precip=0,
         temp=ifelse(treatment=='warming'|treatment=='warming + nutrient addition', 1.5, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment=='warming', 1, ifelse(treatment=='nutrient addition', 3, 4))))%>%
  mutate(resource_mani=ifelse(treatment=='warming', 0, 1))%>%
  mutate(public=0)%>%
  unique()

face<-read.delim("GVN_FACE.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0, 
         p=0, 
         k=0,
         CO2=ifelse(treatment=='A', 0, 160),
         precip=0, 
         temp=0,
         mow_clip=1, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='E', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

nde<-read.csv("IMGERS_NDE.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=="N1M0"|treatment=="N1M1", 1,ifelse(treatment=="N2M0"|treatment=="N2M1", 2, ifelse(treatment=="N3M0"|treatment=="N3M1", 3, ifelse(treatment=="N4M0"|treatment=="N4M1",5, ifelse(treatment=="N5M0"|treatment=="N5M1", 10, ifelse(treatment=="N6M0"|treatment=="N6M1",15, ifelse(treatment=="N7M0"|treatment=="N7M1", 20, ifelse(treatment=="N8M0"|treatment=="N8M1",50,0)))))))), 
         p=0, 
         k=0, 
         CO2=0,
         precip=0, 
         temp=0, 
         mow_clip=ifelse(treatment=="N0M1"|treatment=="N1M1"|treatment=="N2M1"|treatment=="N3M1"|treatment=="N4M1"|treatment=="N5M1"|treatment=="N6M1"|treatment=="N7M1"|treatment=="N8M1", 1,0), 
         burn=0, 
         herb_removal=0, 
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=="N0M0", 0, ifelse(treatment=="N1M0"|treatment=='N0M1'|treatment=="N1M0"|treatment=="N2M0"|treatment=="N3M0"|treatment=="N4M0"|treatment=="N5M0"|treatment=="N6M0"|treatment=="N7M0"|treatment=="N8M0", 1, 2)))%>%
  mutate(resource_mani=ifelse(treatment=="N0M1", 0, 1))%>%
  mutate(public=0)%>%
  unique()

yu<-read.delim("IMGERS_Yu.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
        nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
        n=ifelse(treatment=='N2', 5.6, ifelse(treatment=='N3', 11.2, ifelse(treatment=='N4', 22.4, ifelse(treatment=='N5', 39.2, ifelse(treatment=='N6', 56, 0))))),
        p=ifelse(treatment=='N0', 0, 1.55),
        k=ifelse(treatment=='N0', 0, 3.95),
        CO2=0, 
        precip=0, 
        temp=0,
        mow_clip=1, 
        burn=0, 
        herb_removal=0,
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0, 
        pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0', 0, ifelse(treatment=='N1', 2, 3)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

study119<-read.delim("JRN_Study119.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='T', 10, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='T', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

study278<-read.delim("JRN_study278.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='P1N1'|treatment=='P2N1'|treatment=='P3N1'|treatment=='P4N1'|treatment=='P5N1', 10, 0),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='P1N0'|treatment=='P1N1', -80, ifelse(treatment=='P2N0'|treatment=='P2N1', -50, ifelse(treatment=='P4N0'|treatment=='P4N1', 50, ifelse(treatment=='P5N0'|treatment=='P5N1', 80, 0)))),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='P3N0', 0, ifelse(treatment=='P1N0'|treatment=='P2N0'|treatment=='P3N1'|treatment=='P4N0'|treatment=='P5N0', 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

gce<-read.delim("JSP_GCE2.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=1, water=1, other_manipulation=1,
         n=ifelse(treatment=='N'|treatment=='RN'|treatment=='HN'|treatment=='HRN'|treatment=='CN'|treatment=='CRN'|treatment=='CHN'|treatment=='CHRN', 7, 0),
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='C'|treatment=='CN'|treatment=='CR'|treatment=='CRN'|treatment=='CH'|treatment=='CHN'|treatment=='CHR'|treatment=='CHRN', 300, 0),
         precip=ifelse(treatment=='R'|treatment=='RN'|treatment=='HR'|treatment=='HRN'|treatment=='CR'|treatment=='CRN'|treatment=='CHR'|treatment=='CHRN', 50, 0),
         temp=ifelse(treatment=='H'|treatment=='HN'|treatment=='HR'|treatment=='HRN'|treatment=='CH'|treatment=='CHN'|treatment=='CHR'|treatment=='CHRN', 1.5, 0),
         mow_clip=0,
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='amb', 0, ifelse(treatment=='N'|treatment=='R'|treatment=='H'|treatment=='C', 1, ifelse(treatment=='HRN'|treatment=='CRN'|treatment=='CHN'|treatment=='CHR', 3, ifelse(treatment=='CHRN', 4, 2)))))%>%
  mutate(resource_mani=ifelse(treatment=='H', 0, 1))%>%
  mutate(public=0)%>%
  unique()

wapaclip<-read.delim("KAEFS_WaPaClip.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='U CH'|treatment=='U WH'|treatment=='C CH'|treatment=='C WH', -50, ifelse(treatment=='U CD'|treatment=='U WD'|treatment=='C CD'|treatment=='C WD', 50, 0)),
         temp=ifelse(treatment=='U WC'|treatment=='U WH'|treatment=='U WD'|treatment=='C WC'|treatment=='C WH'|treatment=='C WD', 3, 0),
         mow_clip=ifelse(treatment=='C CC'|treatment=='C CH'|treatment=='C CD'|treatment=='C WC'|treatment=='C WH'|treatment=='C WD', 1, 0),
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='U CC', 0, ifelse(treatment=='U CH', 1, ifelse(treatment=='U CD', 1, ifelse(treatment=='U WC', 1, ifelse(treatment=='C CC', 1, ifelse(treatment=='C WH', 3, ifelse(treatment=='C WD', 3, 2))))))))%>%
  mutate(resource_mani=ifelse(treatment=='U WC', 0, ifelse(treatment=='C CC', 0, ifelse(treatment=='C WC', 0, 1))))%>%
  mutate(public=0)%>%
  unique()

t7<-read.delim("KBS_T7.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='T0F1'|treatment=='T1F1', 12.3, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=0, 
         other_trt=ifelse(treatment=='T1F0'|treatment=='T1F1', 'tilled', 0),
         trt_details=0,
         successional=1, plant_mani=0, pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='T0F0', 0, ifelse(treatment=='T1F1', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='T1F0', 0, 1))%>%
  mutate(public=1)%>%
  unique()

bffert<-read.delim("KLU_BFFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N1F0'|treatment=='N1F1', 17.5, 0),
         p=ifelse(treatment=='N1F0'|treatment=='N1F1', 5, 0),
         k=ifelse(treatment=='N1F0'|treatment=='N1F1', 1.5, 0),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0,
         herb_removal=ifelse(treatment=='N0F1'|treatment=='N1F1', 1, 0),
         trt_details=0,
         other_trt=0, 
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0F0', 0, ifelse(treatment=='N1F1', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='N0F0', 0, ifelse(treatment=='N1F0',3, ifelse(treatment=='N1F1', 4, 1))))%>%
  mutate(public=0)%>%
  unique()

kgfert<-read.delim("KLU_KGFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N1B0'|treatment=='N1B1', 17.5, 0),
         p=ifelse(treatment=='N1B0'|treatment=='N1B1', 5.8, 0),
         k=ifelse(treatment=='N1B0'|treatment=='N1B1', 5.8, 0),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         other_trt=ifelse(treatment=='N0B1'|treatment=='N1B1', "fungicide added", 0), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0B0', 0, ifelse(treatment=='N1B1', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='N0B0', 0, ifelse(treatment=='N1B0',3, ifelse(treatment=='N1B1', 4, 1))))%>%
  mutate(public=0)%>%
  unique()

bgp<-read.delim("KNZ_BGP.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='u_u_p'|treatment=='u_u_c'|treatment=='u_m_p'|treatment=='u_m_c'|treatment=='b_u_p'|treatment=='b_u_c'|treatment=='b_m_p'|treatment=='b_m_c', 0, 10),
         p=ifelse(treatment=='u_u_n'|treatment=='u_u_c'|treatment=='u_m_n'|treatment=='u_m_c'|treatment=='b_u_n'|treatment=='b_u_c'|treatment=='b_m_n'|treatment=='b_m_c', 0, 1),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment=='u_u_n'|treatment=='u_u_p'|treatment=='u_u_c'|treatment=='u_u_b'|treatment=='b_u_n'|treatment=='b_u_p'|treatment=='b_u_c'|treatment=='b_u_b', 0, 1),
         burn=ifelse(treatment=='u_u_n'|treatment=='u_u_p'|treatment=='u_u_c'|treatment=='u_u_b'|treatment=='u_m_n'|treatment=='u_m_p'|treatment=='u_m_c'|treatment=='u_m_b', 0, 1),
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='u_u_c', 0, ifelse(treatment=='u_u_n'|treatment=='u_u_p'|treatment=='u_m_c'|treatment=='b_u_c', 1, ifelse(treatment=='b_u_b'|treatment=='b_m_n'|treatment=='b_m_p'|treatment=='u_m_b', 3, ifelse(treatment=='b_u_n'|treatment=='b_u_p'|treatment=='u_u_b'|treatment=='u_m_n'|treatment=='u_m_p'|treatment=='b_m_c', 2, 4)))))%>%
  mutate(resource_mani=ifelse(treatment=='u_m_c'|treatment=='b_u_c'|treatment=='b_m_c', 0, 1))%>%
  mutate(public=1)%>%
  unique()

irg<-read.delim("KNZ_IRG.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='i', 30, 0),
         temp=0, 
         mow_clip=0, 
         burn=1, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0,
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='i', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

gfp<-read.csv("KNZ_KNP_GFP.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
        n=0, 
        p=0, 
        k=0, 
        CO2=0,
        precip=ifelse(treatment=='Rainout_Grazed'|treatment=="Rainout_Ungrazed", -50, 0),
        temp=0, 
        mow_clip=ifelse(treatment=="Rainout_Ungrazed"|treatment=="Open_Ungrazed",0,1), 
        burn=ifelse(community_type=="20B", 0, 1), 
        herb_removal=ifelse(site_code=="KNP", 1, 0),
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0, 
        pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Open_Ungrazed', 0, ifelse(treatment=="Rainout_Grazed",2,1)))%>%
  mutate(resource_mani=ifelse(treatment=="Open_Grazed",0,1))%>%
  mutate(public=0)%>%
  unique()

pplots<-read.csv("KNZ_PPLOTS.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N1P0'|treatment=='N1P1'|treatment=='N1P2'|treatment=='N1P3', 0, 10),
         p=ifelse(treatment=='N1P1'|treatment=='N2P1', 2.5, ifelse(treatment=='N1P2'|treatment=='N2P2', 5, ifelse(treatment=='N1P3'|treatment=='N2P3', 10, 0))),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N1P0', 0, ifelse(treatment=='N1P1'|treatment=='N1P2'|treatment=='N1P3'|treatment=='N2P0', 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

ramps<-read.csv("KNZ_Ramps.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=0,
         temp=ifelse(treatment=='ambient_heated'|treatment=="delayed_heated", 1, 0),
         mow_clip=0, 
         burn=1, 
         herb_removal=0,
         other_trt=ifelse(treatment=='delayed_control'|treatment=="delayed_heated",'increased precip vari', 'ambient'), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='ambient_control', 0, ifelse(treatment=='delayed_heated', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='ambient_heated', 0, 1))%>%
  mutate(public=1)%>%
  unique()

rhps<-read.delim("KNZ_RHPs.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N'|treatment=='stone+N', 5, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=1, 
         pulse=0,
         other_trt=ifelse(treatment=='stone'|treatment=='stone+N', 'shallow soil', 0),
         trt_details=0,
         successional=1, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment=='stone+N', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='stone', 0, 1))%>%
  mutate(public=1)%>%
  unique()

e6<-read.delim("KUFS_E6.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N0P0S0'|treatment=='N0P8S0', 0, ifelse(treatment=='N4P0S0'|treatment=='N4P8S0', 4, ifelse(treatment=='N8P0S0'|treatment=='N8P8S0', 8, 16))),
         p=ifelse(treatment=='N0P0S0'|treatment=='N4P0S0'|treatment=='N8P0S0', 0, 8),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0P0S0', 0, ifelse(treatment=='N4P0S0'|treatment=='N8P0S0'|treatment=='N16P0S0'|treatment=='N0P8S0', 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

clip<-read.delim("LATNJA_CLIP.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N'|treatment=='TN', 5, 0),
         p=ifelse(treatment=='N'|treatment=='TN', 5, 0), 
         k=0, 
         CO2=0, 
         precip=0,
         temp=ifelse(treatment=='T'|treatment=='TN', 2, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CONTROL', 0, ifelse(treatment=='TN', 3, ifelse(treatment=='T',1,2))))%>%
  mutate(resource_mani=ifelse(treatment=='T', 0, 1))%>%
  mutate(public=0)%>%
  unique()

pme<-read.csv("LEFT_PME.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0,
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=="winwet",50, ifelse(treatment=="winwet_sumwet", 100, 0)),
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=1,
         other_trt=ifelse(treatment=='control', 'ambient precip', ifelse(treatment=='winwet', 'increase winter precip', ifelse(treatment=='winwet_sumdry','increase winter precip, decrease summer precip', ifelse(treatment=='winwet_sumwet', 'increase winter and summer precip', 'decrease winter precip increase summer precip')))), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, 1))%>%#We are considering this to be 1 manipulation, even when they manipulated precip in the winter and summer .
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

herbwood<-read.delim("LG_HerbWood.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='F'|treatment=='FW', 2.4, 0),
         p=ifelse(treatment=='F'|treatment=='FW', 0.66, 0),
         k=ifelse(treatment=='F'|treatment=='FW', 1.5, 0),
         CO2=0,
         precip=ifelse(treatment=='W'|treatment=='FW', 18, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='W', 1, ifelse(treatment=='F', 3, 4))))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

fireplots<-read.delim("MAERC_fireplots.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='snpu'|treatment=='snuu'|treatment=='unpu'|treatment=='unuu'|treatment=='wnpg'|treatment=='wnpu'|treatment=='wnug'|treatment=='wnuu', 5, 0),
         p=ifelse(treatment=='snpu'|treatment=='supu'|treatment=='unpu'|treatment=='uupu'|treatment=='wnpg'|treatment=='wnpu'|treatment=='wupg'|treatment=='wupu', 2, 0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0,
         burn=ifelse(treatment=='uuuu'|treatment=='uupu'|treatment=='unpu'|treatment=='unuu', 0, 1),
         herb_removal=ifelse(treatment=='wnpg'|treatment=='wnug'|treatment=='wupg'|treatment=='wuug', 0, 1), 
         other_trt=0,
         trt_details=ifelse(treatment=='snpu'|treatment=='snuu'|treatment=='supu'|treatment=='suuu', 'summer burn', ifelse(treatment=='wnpg'|treatment=='wnpu'|treatment=='wnug'|treatment=='wnuu'|treatment=='wupg'|treatment=='wupu'|treatment=='wuug'|treatment=='wuuu', 'winter burn', 0)),
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=="wuug", 0, ifelse(treatment=="wnug"|treatment=="wupg"|treatment=="wuuu", 1, ifelse(treatment=="suuu"|treatment=='uuuu'|treatment=='wnpg'|treatment=='wnuu'|treatment=="wupu", 2, ifelse(treatment=='snuu'|treatment=='supu'|treatment=='unuu'|treatment=='uupu'|treatment=='wnpu',3, 4)))))%>%
  mutate(resource_mani=ifelse(treatment=='uuuu'|treatment=='wuuu'|treatment=='suuu', 0, 1))%>%
  mutate(public=0)%>%
  unique()

mwatfer<-read.csv("MNR_watfer.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='F'|treatment=='FW', 10, 0),
         p=ifelse(treatment=='F'|treatment=='FW', 10, 0),
         k=ifelse(treatment=='F'|treatment=='FW', 10, 0),
         CO2=0,
         precip=ifelse(treatment=='W'|treatment=='FW', 18, 0),
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='W', 1, ifelse(treatment=='F', 3, 4))))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

wet<-read.delim("NANT_wet.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='1N0P'|treatment=='1N1P', 67.2, 0),
         p=ifelse(treatment=='0N0P'|treatment=='1N0P', 0, 33.6),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='0N0P', 0, ifelse(treatment=='1N1P', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

gb<-read.delim("NGBER_gb.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=ifelse(treatment=='AMBIENT', 'ambient rainfall', ifelse(treatment=='CURRENT', 'current pattern', ifelse(treatment=='SPRING', 'spring addition', 'winter addition'))), 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='AMBIENT', 0, 1))%>%
  mutate(resource_mani=ifelse(treatment=='CURRENT', 0, 1))%>%
  mutate(public=0)%>%
  unique()

herbdiv<-read.csv("NIN_herbdiv.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='1NF'|treatment=='2NF'|treatment=='3NF'|treatment=='4NF'|treatment=='5NF', 0, 12),
         p=ifelse(treatment=='1NF'|treatment=='2NF'|treatment=='3NF'|treatment=='4NF'|treatment=='5NF', 0, 3.3),
         k=ifelse(treatment=='1NF'|treatment=='2NF'|treatment=='3NF'|treatment=='4NF'|treatment=='5NF', 0, 8),
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=1, 
         herb_removal=ifelse(treatment=='1NF'|treatment=='1F', 0, 1),
         other_trt=0,
         trt_details=ifelse(treatment=='2NF'|treatment=='2F', 'aboveground exclosure', ifelse(treatment=='3NF'|treatment=='3F', 'insecticide', ifelse(treatment=='4NF'|treatment=='4F', 'aboveground exclosure/insecticide', ifelse(treatment=='5NF'|treatment=='5F', 'above/below exclosure/insecticide', 0)))), 
         successional=0, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='1NF', 0, ifelse(treatment=='2NF'|treatment=='3NF', 1, ifelse(treatment=='4NF', 2, ifelse(treatment=='2F'|treatment=='3F', 4, ifelse(treatment=='4F', 5, ifelse(treatment=='1F'|treatment=='5NF', 3, 6)))))))%>%
  mutate(resource_mani=ifelse(treatment=='2NF'|treatment=='3NF'|treatment=='4NF'|treatment=='5NF', 0, 1))%>%
  mutate(public=0)%>%
  unique()

ccd<-read.delim("NTG_CCD.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='CH-'|treatment=='CL-'|treatment=='CN-'|treatment=='WH-'|treatment=='WL-'|treatment=='WN-', -60, ifelse(treatment=='CH+'|treatment=='CL+'|treatment=='CN+'|treatment=='WH+'|treatment=='WL+'|treatment=='WN+', 60, 0)),
         temp=ifelse(site_code=='Alberta'&treatment=='WH-'|site_code=='Alberta'&treatment=='WHA'|site_code=='Alberta'&treatment=='WH+'|site_code=='Alberta'&treatment=='WL-'|site_code=='Alberta'&treatment=='WLA'|site_code=='Alberta'&treatment=='WL+'|site_code=='Alberta'&treatment=='WN-'|site_code=='Alberta'&treatment=='WNA'|site_code=='Alberta'&treatment=='WN+', 2.9,ifelse(site_code=='Manitoba'&treatment=='WH-'|site_code=='Manitoba'&treatment=='WHA'|site_code=='Manitoba'&treatment=='WH+'|site_code=='Manitoba'&treatment=='WL-'|site_code=='Manitoba'&treatment=='WLA'|site_code=='Manitoba'&treatment=='WL+'|site_code=='Manitoba'&treatment=='WN-'|site_code=='Manitoba'&treatment=='WNA'|site_code=='Manitoba'&treatment=='WN+', 1.4,ifelse(site_code=='Saskatchewan'&treatment=='WH-'|site_code=='Saskatchewan'&treatment=='WHA'|site_code=='Saskatchewan'&treatment=='WH+'|site_code=='Saskatchewan'&treatment=='WL-'|site_code=='Saskatchewan'&treatment=='WLA'|site_code=='Saskatchewan'&treatment=='WL+'|site_code=='Saskatchewan'&treatment=='WN-'|site_code=='Saskatchewan'&treatment=='WNA'|site_code=='Saskatchewan'&treatment=='WN+', 1.3,0))),
         mow_clip=ifelse(treatment=='CN-'|treatment=='CNA'|treatment=='CN+'|treatment=='WN-'|treatment=='WN+'|treatment=='WNA', 0, 1),
         burn=0, 
         herb_removal=0, 
         pulse=0,
         other_trt=0,
         trt_details=ifelse(treatment=='CH-'|treatment=='CHA'|treatment=='WH-'|treatment=='CH+'|treatment=='WHA'|treatment=='WH+', 'high intensity defoliation', ifelse(treatment=='CL-'|treatment=='CLA'|treatment=='CL+'|treatment=='WL-'|treatment=='WLA'|treatment=='WL+', 'low intensity defoliation', 0)),
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='CHA'|treatment=='CLA'|treatment=='CN-'|treatment=='CN+', 1,ifelse(treatment=='CNA', 0,ifelse(treatment=='CH-'|treatment=='CH+'|treatment=='CL-'|treatment=='CL+'|treatment=='WHA'|treatment=='WLA'|treatment=='WNA'|treatment=='WN-'|treatment=='WN+', 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='CHA'|treatment=='CLA'|treatment=='WHA'|treatment=='WLA'|treatment=='WNA', 0, 1))%>%
  mutate(public=0)%>%
  unique()

nfert<-read.delim("NWT_246NFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='x', 0, ifelse(treatment=='low', 2, ifelse(treatment=='med', 4, 6))),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='x', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

bowman<-read.delim("NWT_bowman.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N'&calendar_year<=1991, 25, ifelse(treatment=='NP'&calendar_year<=1991, 25, ifelse(treatment=='Control'|treatment=='P', 0, 10))),
         p=ifelse(treatment=='P'&calendar_year<=1991, 25, ifelse(treatment=='NP'&calendar_year<=1991, 25, ifelse(treatment=='Control'|treatment=='N', 0, 10))),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0,
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Control', 0, ifelse(treatment=='NP', 2, 1)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

snow<-read.delim("NWT_snow.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=1,
         n=ifelse(treatment=='XNX'&calendar_year<2011, 28, ifelse(treatment=='XNW'&calendar_year<2011, 28, ifelse(treatment=='PNX'&calendar_year<2011, 28, ifelse(treatment=='PNW'&calendar_year<2011, 28, ifelse(treatment=='XNX'&calendar_year>2010, 10, ifelse(treatment=='XNW'&calendar_year>=2011, 10, ifelse(treatment=='PNX'&calendar_year>=2011, 10, ifelse(treatment=='PNW'&calendar_year>=2011, 10, 0)))))))),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='XXX'|treatment=='XXW'|treatment=='XNX'|treatment=='XNW', 0, 116),
         temp=ifelse(treatment=='XXW'|treatment=='XNW'|treatment=='PXW'|treatment=='PNW', 1, 0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='XXX', 0, ifelse(treatment=='XXW'|treatment=='XNX'|treatment=='PXX', 1, ifelse(treatment=='XNW'|treatment=='PXW'|treatment=='PNX', 2, 3))))%>%
  mutate(resource_mani=ifelse(treatment=='XXW', 0, 1))%>%
  mutate(public=1)%>%
  unique()

oface<-read.delim("ORNL_FACE.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0, 
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='elevated', 170, 0),
         precip=0, 
         temp=0, 
         mow_clip=0, 
         burn=0, 
         herb_removal=0, 
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='elevated', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

tide<-read.delim("PIE_Tide.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='Enriched', 37.5, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='Enriched', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

interaction<-read.delim("RIO_interaction.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N1W1'|treatment=='N1W2'|treatment=='N1W0', 5, 0),
         p=0, 
         k=0,
         CO2=0,
         precip=ifelse(treatment=='N0W0'|treatment=='N1W0'|treatment=='control', 0, 27),
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0,
         trt_details=ifelse(treatment=='N0W1'|treatment=='N1W1', 'small precip pulse', ifelse(treatment=='N0W2'|treatment=='N1W2', 'large precip pulse', 0)),
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='control', 0, ifelse(treatment=='N1W0'|treatment=='N0W1'|treatment=='N0W2', 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

lucero<-read.csv("SCL_Lucero.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N1', 20, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N1', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

ter<-read.csv("SCL_TER.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='OF'|treatment=='CF', 20, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment=='CO'|treatment=='CF', 1, 0), 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='OO', 0, ifelse(treatment=="CF",2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=="CO", 0, 1))%>%
  mutate(public=0)%>%
  unique()

cxn<-read.csv("SERC_CXN.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
        nutrients=1, light=0, carbon=1, water=0, other_manipulation=0,
        n=ifelse(treatment=='t2'|treatment=='t4', 25, 0),
        p=0, 
        k=0, 
        CO2=ifelse(treatment=='t3'|treatment=='t4', 340,0),
        precip=0, 
        temp=0,
        mow_clip=0, 
        burn=0, 
        herb_removal=0,
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0, 
        pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='t1', 0, ifelse(treatment=='t4',2,1)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

tmece<-read.csv("SERC_TMECE.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
        n=0,
        p=0, 
        k=0, 
        CO2=ifelse(treatment=='E', 340,0), 
        precip=0, 
        temp=0,
        mow_clip=0, 
        burn=0, 
        herb_removal=0,
        other_trt=0, 
        trt_details=0,
        successional=0, 
        plant_mani=0, 
        pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='A', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

snfert<-read.delim("SEV_NFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='F', 10, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='F', 1, 0))%>%
  mutate(resource_mani=1)%>%
  mutate(public=1)%>%
  unique()

wenndex<-read.delim("SEV_WENNDEx.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=1,
         n=ifelse(treatment=='C'|treatment=='P'|treatment=='T'|treatment=='TP', 0, 2),
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='C'|treatment=='N'|treatment=='T'|treatment=='TN', 0, 50),
         temp=ifelse(treatment=='C'|treatment=='N'|treatment=='P'|treatment=='PN', 0, 1),
         mow_clip=0, 
         burn=1, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='N'|treatment=='P'|treatment=='T', 1, ifelse(treatment=='TPN', 3, 2))))%>%
  mutate(resource_mani=ifelse(treatment=='T', 0, 1))%>%
  mutate(public=1)%>%
  unique()

grazeprecip<-read.csv("SFREC_GrazePrecip.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0,
         p=0, 
         k=0, 
         CO2=0, 
         precip=ifelse(treatment=="C", 0, ifelse(treatment=='W',80,-80)), 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

uk<-read.delim("SKY_UK.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=0,
         precip=ifelse(treatment=='C'|treatment=='H', 0, 30),
         temp=ifelse(treatment=='C'|treatment=='P', 0, 3),
         mow_clip=1, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=1, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='HP', 2, 1)))%>%
  mutate(resource_mani=ifelse(treatment=='H', 0, 1))%>%
  mutate(public=0)%>%
  unique()

nitrogen<-read.csv("SR_Nitrogen.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='1_NITROGEN'|treatment=='0_NITROGEN', 4, 0),
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=ifelse(treatment=='1_CONTROL'|treatment=='1_NITROGEN',1,0), 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='0_CONTROL', 0, ifelse(treatment=='1_NITROGEN',2,1)))%>%
  mutate(resource_mani=ifelse(treatment=='1_CONTROL',0,1))%>%
  mutate(public=0)%>%
  unique()

water<-read.csv("SR_Water.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment, community_type)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0,
         p=ifelse(treatment=='0_WATER_0'|treatment=='0_WATER_1'|treatment=='1_WATER_0'|treatment=='1_WATER_1', 34.1,0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=ifelse(treatment=='0_CONTROL_0'|treatment=='1_CONTROL_0'|treatment=="0_WATER_0"|treatment=='1_WATER_0', 1,0),
         other_trt=0, 
         trt_details=0,
         successional=1, 
         plant_mani=ifelse(treatment=='1_CONTROL_0'|treatment=='1_CONTROL_1'|treatment=="1_WATER_0"|treatment=='1_WATER_1',1,0), pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='0_CONTROL_1', 0, ifelse(treatment=='1_WATER_0',3,ifelse(treatment=='1_CONTROL_1'|treatment=='0_CONTROL_0'|treatment=='0_WATER_1',1,2))))%>%
  mutate(resource_mani=ifelse(treatment=='1_CONTROL_1'|treatment=='0_CONTROL_0'|treatment=="1_CONTROL_0",0,1))%>%
  mutate(public=0)%>%
  unique()

gane<-read.delim("SVA_GANE.txt")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0, 
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='C'|treatment=='P', 0, ifelse(treatment=='LN'|treatment=='LNP', 0.5, 5)),
         p=ifelse(treatment=='P'|treatment=='LNP'|treatment=='HNP', 1, 0),
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=1,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='C', 0, ifelse(treatment=='LN'|treatment=='HN'|treatment=='P', 1, 2)))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

tface<-read.csv("TAS_FACE.csv")%>%
 select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0,
         nutrients=0, light=0, carbon=1, water=0, other_manipulation=1,
         n=0, 
         p=0, 
         k=0, 
         CO2=ifelse(treatment=='UnwarmedFACE'|treatment=="WarmedFACE", 170,0), 
         precip=0, 
         temp=ifelse(treatment=='WarmedControl'|treatment=="WarmedFACE",2,0),
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='UnwarmedControl', 0, ifelse(treatment=='WarmedFACE',2,1)))%>%
  mutate(resource_mani=ifelse(treatment=='WarmedControl',0,1))%>%
  mutate(public=0)%>%
  unique()

lovegrass<-read.csv("TRA_Lovegrass.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0,
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='gcc'|treatment=='ghc'|treatment=='gsc'|treatment=="ncc"|treatment=='nhc'|treatment=='nsc', 0, 0.432), 
         p=ifelse(treatment=='gcc'|treatment=='ghc'|treatment=='gsc'|treatment=="ncc"|treatment=='nhc'|treatment=='nsc', 0, 0.022), 
         k=ifelse(treatment=='gcc'|treatment=='ghc'|treatment=='gsc'|treatment=="ncc"|treatment=='nhc'|treatment=='nsc', 0, 0.082), 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=ifelse(treatment=='gsc'|treatment=='gsn'|treatment=='nsc'|treatment=="nsn",1, 0), 
         burn=0, 
         herb_removal=ifelse(treatment=='ncc'|treatment=='ncn'|treatment=='nhc'|treatment=="nhn"|treatment=="nsc"|treatment=='nsn',1,0),
         other_trt=0,
         trt_details=0, 
         successional=0, 
         plant_mani=ifelse(treatment=='ghc'|treatment=='ghn'|treatment=='nhc'|treatment=='nhn',1,0), 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='gcc',0,ifelse(treatment=='ghc'|treatment=='gsc'|treatment=='ncc',1, ifelse(treatment=='nhc'|treatment=='nsc',2,ifelse(treatment=="gcn",3,ifelse(treatment=='nhn'|treatment=='nsn',5,4))))))%>%
  mutate(resource_mani=ifelse(treatment=='ghc'|treatment=='gsc'|treatment=="ncc"|treatment=='nhc'|treatment=='nsc', 0, 1))%>%
  mutate(public=0)%>%
  unique()

nitadd<-read.csv("YMN_NitAdd.csv")%>%
  select(site_code, project_name, calendar_year, treatment_year, treatment)%>%
  mutate(community_type=0,
         nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N0',0, ifelse(treatment=='N5', 2.3, ifelse(treatment=='N10',4.7, ifelse(treatment=='N20',9.3, ifelse(treatment=='N40',18.7,37.3))))), 
         p=0, 
         k=0, 
         CO2=0, 
         precip=0, 
         temp=0,
         mow_clip=0, 
         burn=0, 
         herb_removal=0,
         other_trt=0, 
         trt_details=0,
         successional=0, 
         plant_mani=0, 
         pulse=0)%>%
  mutate(plot_mani=ifelse(treatment=='N0', 0, 1))%>%
  mutate(resource_mani=1)%>%
  mutate(public=0)%>%
  unique()

###merge all datasets
combine<-rbind(bffert, bgp, biocon, bowman, ccd, clip, clonal, culardoch, cxn, e001, e002, e6, events, exp1, face, fireplots,gane, gap2, gb, gce, gfp, grazeprecip, herbdiv, herbwood, imagine, interaction, irg, kgfert, lind, lovegrass, lucero, mat2, megarich, mnt, mwatfer, nde, nfert, nitadd, nitphos, nitrogen,nsfc, oface, pennings, pplots,pme, pq, ramps, rhps, rmapc, snfert, snow, study119, study278, t7, ter, tface,tide,tmece, uk, wapaclip, warmnut, water, watering, wenndex, wet, yu)

#kim's
write.csv(combine, 'C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\ExperimentInformation_Feb2016a.csv')

#meghan's
write.csv(combine, "~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_March2016.csv")

# ##checking
# check<-combine%>%
#   select(-calendar_year, -treatment_year)%>%
#   unique()
# write.csv(check, "~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Feb2016_CHECK.csv")
