#notes:
#if clipping, burning are light treatments, then the following experiments need to have light=1: megarich, e001?, wapaclip, bgp

setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\clean datasets - please do not touch\\sp text files')

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/FINAL_SEPT2014/clean datasets - please do not touch/sp text files")

library(tidyr)
library(dplyr)

# nme<-read.delim("Alps_NME.txt")%>%
#   select(site_code, project_name, calendar_year, treatment)%>%
#   mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
#          n=0, p=0, k=0, other_nut=0, lime=0, 
#          soil_carbon=ifelse(treatment=='CG', 15, ifelse(treatment=='CU', 15, 0)), 
#          CO2=0, precip=0, temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0,
#          grazed=ifelse(treatment=='CoG', 1, ifelse(treatment=='CG', 1, 0)),
#          fungicide=0, herb_removal=0, ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
#   unique()

watering<-read.delim("ANG_watering.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='W', 20, ifelse(treatment=='S', 20, 0)), temp=0, precip_vari=0,
         precip_season=ifelse(treatment=='W', 'winter addition', ifelse(treatment=='S', 'spring addition', 0)),
         mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

mat2<-read.delim("ARC_mat2.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NP', 10, 0), p=ifelse(treatment=='NP', 5, 0), k=0, other_nut=0,
         lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

mnt<-read.delim("ARC_mnt.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NP', 10, 0), p=ifelse(treatment=='NP', 5, 0), k=0, other_nut=0,
         lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

clonal<-read.delim("ASGA_Clonal.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='mixed_CO', 0, ifelse(treatment=='non-clonal_CO', 0, 20.1)), p=0, k=0, other_nut=0,
         lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0,
         other_trt=ifelse(treatment=='non-clonal_CO', 'non-clonal species', ifelse(treatment=='mixed_LP', 'large nutrient patches', ifelse(treatment=='non-clonal_LP', 'non-clonal species, large nutrient patches', ifelse(treatment=='mixed_SP', 'small nutrient patches', ifelse(treatment=='non-clonal_SP', 'non-clonal species, small nutrient patches', ifelse(treatment=='non-clonal_UN', 'non-clonal species', 0)))))),
         initial_disturbance=1, plant_mani=1, cessation=0)%>%
  unique()

exp1<-read.delim("ASGA_Exp1.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='2_0_CO', 0, ifelse(treatment=='1_0_CO', 0, ifelse(treatment=='2_1_CO', 0, ifelse(treatment=='1_0_CO', 0, 20.1)))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0,
         other_trt=ifelse(treatment=='2_0_PA', 'nutrient patches', ifelse(treatment=='1_0_PA', 'nutrient patches', ifelse(treatment=='2_1_PA', 'nutrient patches', ifelse(treatment=='1_1_PA', 'nutrient patches', 0)))),
         initial_disturbance=ifelse(treatment=='2_0_CO', 1, ifelse(treatment=='2_0_PA', 1, ifelse(treatment=='2_0_UN', 1, 0))),
         plant_mani=ifelse(treatment=='2_1_CO', 1, ifelse(treatment=='2_1_PA', 1, ifelse(treatment=='2_1_UN', 1, ifelse(treatment=='1_1_CO', 1, ifelse(treatment=='1_1_PA', 1, ifelse(treatment=='1_1_UN', 1, 0)))))),
         cessation=0)%>%
  unique()

#do we want to include what the species treatments are here? do we want to drop anything but the 16 species sowed plots?
lind<-read.delim("BAY_LIND.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='rain_rich1', 8, ifelse(treatment=='rain_rich2', 8, ifelse(treatment=='rain_rich4', 8, ifelse(treatment=='rain_rich8', 8, ifelse(treatment=='rain_rich16', 8, 0))))),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=1, plant_mani=1, cessation=0)%>%
  unique()

events<-read.delim("Bt_EVENT2.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=ifelse(treatment=='CM-N1', 'reduced variability', ifelse(treatment=='D1-N1', 'early drought', ifelse(treatment=='D2-N1', 'late dorught', 0))),
         precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

pq<-read.delim("BUX_PQ.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=ifelse(treatment=='wet', 20, ifelse(treatment=='warm wet', 20, ifelse(treatment=='dry', -20, ifelse(treatment=='warm dry', -20, 0)))),
         temp=ifelse(treatment=='warm', 3, ifelse(treatment=='warm dry', 3, ifelse(treatment=='warm wet', 3, 0))),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

pennings<-read.delim("CAR_Pennings.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='NPK'&calendar_year>1999, 164, ifelse(treatment=='NPK'&calendar_year==1999, 84, 0)),
         p=ifelse(treatment=='NPK'&calendar_year>1999, 82, ifelse(treatment=='NPK'&calendar_year==1999, 42, 0)),
         k=ifelse(treatment=='NPK'&calendar_year>1999, 41, ifelse(treatment=='NPK'&calendar_year==1999, 21, 0)),
         other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

rmapc<-read.delim("CAU_RMAPC.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N', 9, ifelse(treatment=='NP', 9, 0)),
         p=ifelse(treatment=='P', 2.6, ifelse(treatment=='NP', 2.6, 0)),
         k=0, other_nut=0, lime=ifelse(treatment=='Ca', 1, 0), soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='H2O', 20, 0),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

biocon<-read.delim("CDR_biocon.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=1, water=0, other_manipulation=1,
         n=ifelse(treatment=='N_X', 4, ifelse(treatment=='N_C', 4, 0)),
         p=0, k=0, other_nut=0, lime=0,
         soil_carbon=ifelse(treatment=='X_C', 160, ifelse(treatment=='N_C', 160, 0)),
         CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=1, plant_mani=1, cessation=0)%>%
  unique()

#all plots were burned in some years, but not others. how do we want to code this?
e001<-read.delim("CDR_e001_revised.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='2_y_n', 1.02, ifelse(treatment=='3_y_n', 2.04, ifelse(treatment=='4_y_n', 3.40, ifelse(treatment=='5_y_n', 5.44, ifelse(treatment=='6_y_n', 9.52, ifelse(treatment=='7_y_n', 17, ifelse(treatment=='8_y_n', 27.2, 0))))))),
         p=ifelse(treatment=='9_y_n', 0, 4.6), k=ifelse(treatment=='9_y_n', 0, 6.1),
         other_nut=ifelse(treatment=='9_y_n', 0, 1), lime=ifelse(treatment=='9_y_n', 0, 1),
         soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=1, grazed=0, fungicide=0, herb_removal=1,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

mwatfer<-read.csv("MNR_watfer.csv")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='F', 10, ifelse(treatment=='FW', 10, 0)),
         p=ifelse(treatment=='F', 10, ifelse(treatment=='FW', 10, 0)),
         k=ifelse(treatment=='F', 10, ifelse(treatment=='FW', 10, 0)),
         other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='W', 18, ifelse(treatment=='FW', 18, 0)),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#half of the plots were turned into cessation, so it seems that we might want to delete that half from the entire dataset for the first analysis; we need to think about this more
e002<-read.delim("CDR_e002.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='2_f_u_c', 1.02, ifelse(treatment=='2_f_u_n', 1.02, ifelse(treatment=='3_f_u_c', 2.04, ifelse(treatment=='3_f_u_n', 2.04, ifelse(treatment=='4_f_u_c', 3.4, ifelse(treatment=='4_f_u_n', 3.4, ifelse(treatment=='5_f_u_c', 5.44, ifelse(treatment=='5_f_u_n', 5.44, ifelse(treatment=='6_f_u_c', 9.52, ifelse(treatment=='6_f_u_n', 9.52, ifelse(treatment=='7_f_u_c', 17, ifelse(treatment=='7_f_u_n', 17, ifelse(treatment=='8_f_u_c', 27.2, ifelse(treatment=='8_f_u_n', 27.2, 0)))))))))))))),
         p=ifelse(treatment=='9_f_u_c', 0, ifelse(treatment=='9_f_u_n', 0, 4.6)),
         k=ifelse(treatment=='9_f_u_c', 0, ifelse(treatment=='9_f_u_n', 0, 6.1)),
         other_nut=ifelse(treatment=='9_f_u_c', 0, ifelse(treatment=='9_f_u_n', 0, 1)),
         lime=ifelse(treatment=='9_f_u_c', 0, ifelse(treatment=='9_f_u_n', 0, 1)),
         soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=1, plant_mani=0,
         cessation=ifelse(treatment=='1_f_u_c', 1, ifelse(treatment=='2_f_u_c', 1, ifelse(treatment=='3_f_u_c', 1, ifelse(treatment=='4_f_u_c', 1, ifelse(treatment=='5_f_u_c', 1, ifelse(treatment=='6_f_u_c', 1, ifelse(treatment=='7_f_u_c', 1, ifelse(treatment=='8_f_u_c', 1, ifelse(treatment=='9_f_u_c', 1, 0))))))))))%>%
  unique()

#these are the megaliths, so I labeled this as initial disturbance; should we rethink this?
megarich<-read.delim("CEH_Megarich.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=1, water=0, other_manipulation=1,
         n=10, p=2, k=2, other_nut=0, lime=0, soil_carbon=0,
         CO2=ifelse(treatment=='EcAt', 280, ifelse(treatment=='EcEt', 280, 0)),
         precip=0,
         temp=ifelse(treatment=='AcEt', 2.9, ifelse(treatment=='EcEt', 2.9, 0)),
         precip_vari=0, precip_season=0, mow_clip=1, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=1, plant_mani=0, cessation=0)%>%
  unique()

yu<-read.delim("IMGERS_Yu.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N2', 5.6, ifelse(treatment=='N3', 11.2, ifelse(treatment=='N4', 22.4, ifelse(treatment=='N5', 39.2, ifelse(treatment=='N6', 56, 0))))),
         p=ifelse(treatment=='N0', 0, 1.55),
         k=ifelse(treatment=='N0', 0, 3.95),
         other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

interaction<-read.delim("RIO_interaction.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N1W1', 5, ifelse(treatment=='N1W2', 5, ifelse(treatment=='N1W0', 5, 0))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='N0W0', 0, ifelse(treatment=='N1W0', 0, ifelse(treatment=='control', 0, 27))),
         temp=0,
         precip_vari=ifelse(treatment=='N0W1', 'small pulse', ifelse(treatment=='N1W1', 'small pulse', ifelse(treatment=='N0W2', 'large pulse', ifelse(treatment=='N1W2', 'large pulse', 0)))),
         precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

imagine<-read.delim("CLE_imagine.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=1, water=1, other_manipulation=1,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0,
         CO2=ifelse(treatment=='TDCO2', 200, 0),
         precip=ifelse(treatment=='TD', -20, ifelse(treatment=='TDCO2', -20, 0)),
         temp=ifelse(treatment=='C', 0, 3.5),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#one time burn in 1999 for a subset of the plots, do we want to code this by year?
culardoch<-read.delim("CUL_culardoch.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N10', 1, ifelse(treatment=='N10burn', 1, ifelse(treatment=='N10clip', 1, ifelse(treatment=='N10burnclip', 1, ifelse(treatment=='N20', 2, ifelse(treatment=='N20burn', 2, ifelse(treatment=='N20clip', 2, ifelse(treatment=='N20burnclip', 2, ifelse(treatment=='N50', 5, ifelse(treatment=='N50burn', 5, ifelse(treatment=='N50clip', 5, ifelse(treatment=='N50burnclip', 5, 0)))))))))))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0,
         mow_clip=ifelse(treatment=='clip', 1, ifelse(treatment=='burnclip', 1, ifelse(treatment=='N10clip', 1, ifelse(treatment=='N20clip', 1, ifelse(treatment=='N50clip', 1, ifelse(treatment=='N10burnclip', 1, ifelse(treatment=='N20burnclip', 1, ifelse(treatment=='N50burnclip', 1, 0)))))))),
         burn=ifelse(treatment=='N10burn', 1, ifelse(treatment=='N20burn', 1, ifelse(treatment=='N50burn', 1, ifelse(treatment=='burn', 1, ifelse(treatment=='burnclip', 1, ifelse(treatment=='N10burnclip', 1, ifelse(treatment=='N20burnclip', 1, ifelse(treatment=='N50burnclip', 1, 0)))))))),
         grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#coded as initial disturbance because of the pulse nature of the gap openings
gap2<-read.delim("DCGS_gap2.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=1, carbon=0, water=0, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=ifelse(treatment=='_018', '18 ft opening', ifelse(treatment=='_033', '33 ft opening', ifelse(treatment=='_066', '66 ft opening', ifelse(treatment=='_100', '100 ft opening', ifelse(treatment=='_150', '150 ft opening', 0))))), other_trt=0, initial_disturbance=1, plant_mani=0, cessation=0)%>%
  unique()

nsfc<-read.delim("DL_NSFC.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='N', 10, ifelse(treatment=='WN', 10, 0)),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='W', 49.8, ifelse(treatment=='WN', 49.8, 0)),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

warmnut<-read.delim("Finse_WarmNut.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='nutrient addition', 10, ifelse(treatment=='warming + nutrient addition', 10, 0)),
         p=ifelse(treatment=='nutrient addition', 2, ifelse(treatment=='warming + nutrient addition', 2, 0)),
         k=ifelse(treatment=='nutrient addition', 8, ifelse(treatment=='warming + nutrient addition', 8, 0)),
         other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0,
         temp=ifelse(treatment=='warming', 1.5, ifelse(treatment=='warming + nutrient addition', 1.5, 0)),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

face<-read.delim("GVN_FACE.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0,
         CO2=ifelse(treatment=='A', 0, 160),
         precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

study119<-read.delim("JRN_Study119.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='T', 10, 0),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

study278<-read.delim("JRN_study278.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='P1N1', 10, ifelse(treatment=='P2N1', 10, ifelse(treatment=='P3N1', 10, ifelse(treatment=='P4N1', 10, ifelse(treatment=='P5N1', 10, 0))))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='P1N0', -80, ifelse(treatment=='P1N1', -80, ifelse(treatment=='P2N0', -50, ifelse(treatment=='P2N1', -50, ifelse(treatment=='P4N0', 50, ifelse(treatment=='P4N1', 50, ifelse(treatment=='P5N0', 80, ifelse(treatment=='P5N1', 80, 0)))))))),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

gce<-read.delim("JSP_GCE2.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=1, water=1, other_manipulation=1,
         n=ifelse(treatment=='N', 7, ifelse(treatment=='RN', 7, ifelse(treatment=='HN', 7, ifelse(treatment=='HRN', 7, ifelse(treatment=='CN', 7, ifelse(treatment=='CRN', 7, ifelse(treatment=='CHN', 7, ifelse(treatment=='CHRN', 7, 0)))))))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0,
         CO2=ifelse(treatment=='C', 300, ifelse(treatment=='CN', 300, ifelse(treatment=='CR', 300, ifelse(treatment=='CRN', 300, ifelse(treatment=='CH', 300, ifelse(treatment=='CHN', 300, ifelse(treatment=='CHR', 300, ifelse(treatment=='CHRN', 300, 0)))))))),
         precip=ifelse(treatment=='R', 50, ifelse(treatment=='RN', 50, ifelse(treatment=='HR', 50, ifelse(treatment=='HRN', 50, ifelse(treatment=='CR', 50, ifelse(treatment=='CRN', 50, ifelse(treatment=='CHR', 50, ifelse(treatment=='CHRN', 50, 0)))))))),
         temp=ifelse(treatment=='H', 1.5, ifelse(treatment=='HN', 1.5, ifelse(treatment=='HR', 1.5, ifelse(treatment=='HRN', 1.5, ifelse(treatment=='CH', 1.5, ifelse(treatment=='CHN', 1.5, ifelse(treatment=='CHR', 1.5, ifelse(treatment=='CHRN', 1.5, 0)))))))),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#if we consider clipping a light treatment, this would be one
#I don't understand why we dropped part of this before, it seems the clipping is just a treatment
wapaclip<-read.delim("KAEFS_WaPaClip.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='U CH', -50, ifelse(treatment=='U CD', 50, ifelse(treatment=='U WH', -50, ifelse(treatment=='U WD', 50, ifelse(treatment=='C CH', -50, ifelse(treatment=='C CD', 50, ifelse(treatment=='C WH', -50, ifelse(treatment=='C WD', 50, 0)))))))),
         temp=ifelse(treatment=='U WC', 3, ifelse(treatment=='U WH', 3, ifelse(treatment=='U WD', 3, ifelse(treatment=='C WC', 3, ifelse(treatment=='C WH', 3, ifelse(treatment=='C WD', 3, 0)))))),
         precip_vari=0, precip_season=0,
         mow_clip=ifelse(treatment=='C CC', 1, ifelse(treatment=='C CH', 1, ifelse(treatment=='C CD', 1, ifelse(treatment=='C WC', 1, ifelse(treatment=='C WH', 1, ifelse(treatment=='C WD', 1, 0)))))),
         burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#all plots started to be burned in 1997, so we might want to include that
t7<-read.delim("KBS_T7.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='T0F1', 12.3, ifelse(treatment=='T1F1', 12.3, 0)),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0,
         other_trt=ifelse(treatment=='T1F0', 'tilled', ifelse(treatment=='T1F1', 'tilled', 0)),
         initial_disturbance=1, plant_mani=0, cessation=0)%>%
  unique()

bffert<-read.delim("KLU_BFFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N1F0', 17.5, ifelse(treatment=='N1F1', 17.5, 0)),
         p=ifelse(treatment=='N1F0', 5, ifelse(treatment=='N1F1', 5, 0)),
         k=ifelse(treatment=='N1F0', 1.5, ifelse(treatment=='N1F1', 1.5, 0)),
         other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0,
         herb_removal=ifelse(treatment=='N0F1', 1, ifelse(treatment=='N1F1', 1, 0)),
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

kgfert<-read.delim("KLU_KGFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N1B0', 17.5, ifelse(treatment=='N1B1', 17.5, 0)),
         p=ifelse(treatment=='N1B0', 5.8, ifelse(treatment=='N1B1', 5.8, 0)),
         k=ifelse(treatment=='N1B0', 5.8, ifelse(treatment=='N1B1', 5.8, 0)),
         other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0,
         fungicide=ifelse(treatment=='N0B1', 1, ifelse(treatment=='N1B1', 1, 0)),
         herb_removal=0, ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

bgp<-read.delim("KNZ_BGP.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='u_u_p', 0, ifelse(treatment=='u_u_c', 0, ifelse(treatment=='u_m_p', 0, ifelse(treatment=='u_m_c', 0, ifelse(treatment=='b_u_p', 0, ifelse(treatment=='b_u_c', 0, ifelse(treatment=='b_m_p', 0, ifelse(treatment=='b_m_c', 0, 10)))))))),
         p=ifelse(treatment=='u_u_n', 0, ifelse(treatment=='u_u_c', 0, ifelse(treatment=='u_m_n', 0, ifelse(treatment=='u_m_c', 0, ifelse(treatment=='b_u_n', 0, ifelse(treatment=='b_u_c', 0, ifelse(treatment=='b_m_n', 0, ifelse(treatment=='b_m_c', 0, 10)))))))),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0,
         mow_clip=ifelse(treatment=='u_u_n', 0, ifelse(treatment=='u_u_p', 0, ifelse(treatment=='u_u_c', 0, ifelse(treatment=='u_u_b', 0, ifelse(treatment=='b_u_n', 0, ifelse(treatment=='b_u_p', 0, ifelse(treatment=='b_u_c', 0, ifelse(treatment=='b_u_b', 0, 1)))))))),
         burn=ifelse(treatment=='u_u_n', 0, ifelse(treatment=='u_u_p', 0, ifelse(treatment=='u_u_c', 0, ifelse(treatment=='u_u_b', 0, ifelse(treatment=='u_m_n', 0, ifelse(treatment=='u_m_p', 0, ifelse(treatment=='u_m_c', 0, ifelse(treatment=='u_m_b', 0, 1)))))))),
         grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#annual burn, do we include this?
irg<-read.delim("KNZ_IRG.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='i', 30, 0),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#2 year burn, do we include this?
pplots<-read.delim("KNZ_PPLOTS.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N0P0', 0, ifelse(treatment=='N0P1', 0, ifelse(treatment=='N0P2', 0, ifelse(treatment=='N0P3', 0, 10)))),
         p=ifelse(treatment=='N0P1', 2.5, ifelse(treatment=='N0P2', 5, ifelse(treatment=='N0P3', 10, ifelse(treatment=='N1P1', 2.5, ifelse(treatment=='N1P2', 5, ifelse(treatment=='N1P3', 10, 0)))))),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#annual burn, do we include this?
ramps<-read.delim("KNZ_Ramps.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='AD', -30, ifelse(treatment=='DD', -30, 0)),
         temp=ifelse(treatment=='AH', 1, ifelse(treatment=='DH', 1, 0)),
         precip_vari=ifelse(treatment=='DC', 'delayed', ifelse(treatment=='DH', 'delayed', ifelse(treatment=='DD', 'delayed', 'ambient'))),
         precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

rhps<-read.delim("KNZ_RHPs.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N', 5, ifelse(treatment=='stone+N', 5, 0)),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0,
         other_trt=ifelse(treatment=='stone', 'shallow soil', ifelse(treatment=='stone+N', 'shallow soil', 0)),
         initial_disturbance=1, plant_mani=0, cessation=0)%>%
  unique()

e6<-read.delim("KUFS_E6.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N0P0S0', 0, ifelse(treatment=='N0P8S0', 0, ifelse(treatment=='N4P0S0', 4, ifelse(treatment=='N4P8S0', 4, ifelse(treatment=='N8P0S0', 8, ifelse(treatment=='N8P8S0', 8, 16)))))),
         p=ifelse(treatment=='N0P0S0', 0, ifelse(treatment=='N4P0S0', 0, ifelse(treatment=='N8P0S0', 0, 8))),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

clip<-read.delim("LATNJA_CLIP.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='N', 5, ifelse(treatment=='TN', 5, 0)),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0,
         temp=ifelse(treatment=='T', 2, ifelse(treatment=='TN', 2, 0)),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

herbwood<-read.delim("LG_HerbWood.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='F', 2.4, ifelse(treatment=='FW', 2.4, 0)),
         p=ifelse(treatment=='F', 0.66, ifelse(treatment=='FW', 0.66, 0)),
         k=ifelse(treatment=='F', 1.5, ifelse(treatment=='FW', 1.5, 0)),
         other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='W', 18, ifelse(treatment=='FW', 18, 0)),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#they split into summer and winter burning, where should we code that? all other experiment's burn treatments are just binary
fireplots<-read.delim("MAERC_fireplots.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=1,
         n=ifelse(treatment=='snpg', 5, ifelse(treatment=='snpu', 5, ifelse(treatment=='snug', 5, ifelse(treatment=='snuu', 5, ifelse(treatment=='unpg', 5, ifelse(treatment=='unpu', 5, ifelse(treatment=='unug', 5, ifelse(treatment=='unuu', 5, ifelse(treatment=='wnpg', 5, ifelse(treatment=='wnpu', 5, ifelse(treatment=='wnug', 5, ifelse(treatment=='wnuu', 5, 0)))))))))))),
         p=ifelse(treatment=='snpg', 5, ifelse(treatment=='snpu', 5, ifelse(treatment=='supg', 5, ifelse(treatment=='supu', 5, ifelse(treatment=='unpg', 5, ifelse(treatment=='unpu', 5, ifelse(treatment=='uupg', 5, ifelse(treatment=='uupu', 5, ifelse(treatment=='wnpg', 5, ifelse(treatment=='wnpu', 5, ifelse(treatment=='wupg', 5, ifelse(treatment=='wupu', 5, 0)))))))))))),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0,
         burn=ifelse(treatment=='uuuu', 0, ifelse(treatment=='uupu', 0, ifelse(treatment=='unpu', 0, ifelse(treatment=='uuug', 0, ifelse(treatment=='unug', 0, ifelse(treatment=='uupg', 0, ifelse(treatment=='unpg', 0, ifelse(treatment=='unuu', 0, 1)))))))),
         grazed=ifelse(treatment=='snpg', 1, ifelse(treatment=='snug', 1, ifelse(treatment=='supg', 1, ifelse(treatment=='suug', 1, ifelse(treatment=='unpg', 1, ifelse(treatment=='unug', 1, ifelse(treatment=='uupg', 1, ifelse(treatment=='uuug', 1, ifelse(treatment=='wnpg', 1, ifelse(treatment=='wnug', 1, ifelse(treatment=='wupg', 1, ifelse(treatment=='wuug', 1, 0)))))))))))),
         fungicide=0, herb_removal=0, ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

wet<-read.delim("NANT_wet.txt")%>%
  select(site_code, project_name, community_type, calendar_year, treatment)%>%
  mutate(nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='1N0P', 67.2, ifelse(treatment=='1N1P', 67.2, 0)),
         p=ifelse(treatment=='0N0P', 0, ifelse(treatment=='1N0P', 0, 33.6)),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

gb<-read.delim("NGBER_gb.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0,
         precip_season=ifelse(treatment=='AMBIENT', 0, ifelse(treatment=='CURRENT', 'current pattern', ifelse(treatment=='SPRING', 'spring addition', 'winter addition'))),
         mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#I had been thinking plot_mani would just be calculated by making all treatment columns binary, and then summing across them; but this one has so many additional treatments, so I'm not sure if I should be going back to add in a plot-mani column and manually putting in what we think plot mani is for each treatment. There must be a better way.
# herbdiv<-read.delim("NIN_herbdiv.txt")%>%
#   select(site_code, project_name, calendar_year, treatment)%>%
#   mutate(community_type=0, nutrients=0, light=0, carbon=0, water=0, other_manipulation=0,
#          n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
#          precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
#          ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
#   unique()

# this should really be split into three sites, because they are provinces across canada
ccd<-read.delim("NTG_CCD.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=0, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

nfert<-read.delim("NWT_246NFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='x', 0, ifelse(treatment=='low', 2, ifelse(treatment=='med', 4, 6))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

bowman<-read.delim("NWT_bowman.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='N'&calendar_year==1990, 25, ifelse(treatment=='NP'&calendar_year==1990, 25, ifelse(treatment=='N'&calendar_year==1991, 0, ifelse(treatment=='NP'&calendar_year==1991, 0, ifelse(treatment=='Control', 0, ifelse(treatment=='P', 0, 10)))))),
         p=ifelse(treatment=='P'&calendar_year==1990, 25, ifelse(treatment=='NP'&calendar_year==1990, 25, ifelse(treatment=='P'&calendar_year==1991, 0, ifelse(treatment=='NP'&calendar_year==1991, 0, ifelse(treatment=='Control', 0, ifelse(treatment=='N', 0, 10)))))),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

snow<-read.delim("NWT_snow.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=1,
         n=ifelse(treatment=='XNX'&calendar_year<2011, 28, ifelse(treatment=='XNW'&calendar_year<2011, 28, ifelse(treatment=='PNX'&calendar_year<2011, 28, ifelse(treatment=='PNW'&calendar_year<2011, 28, ifelse(treatment=='XNX'&calendar_year>2010, 10, ifelse(treatment=='XNW'&calendar_year>=2011, 10, ifelse(treatment=='PNX'&calendar_year>=2011, 10, ifelse(treatment=='PNW'&calendar_year>=2011, 10, 0)))))))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='XXX', 0, ifelse(treatment=='XXW', 0, ifelse(treatment=='XNX', 0, ifelse(treatment=='XNW', 0, 116)))),
         temp=ifelse(treatment=='XXW', 1, ifelse(treatment=='XNW', 1, ifelse(treatment=='PXW', 1, ifelse(treatment=='PNW', 1, 0)))),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#this includes woody and herbaceous species covers
#biomass includes woody and herbaceous species from 2001-2003, but only herbaceous in 2008
oface<-read.delim("ORNL_FACE.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=1, water=0, other_manipulation=0,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0,
         CO2=ifelse(treatment=='elevated', 170, 0),
         precip=0, temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0,
         herb_removal=0, ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#this has two community types, one of which got higher nutrients... we average them before
tide<-read.delim("PIE_Tide.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='Enriched', 37.5, 0),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

snfert<-read.delim("SEV_NFert.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='F', 10, 0),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#wildfire burnt all plots in 2009
wenndex<-read.delim("SEV_WENNDEx.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=1,
         n=ifelse(treatment=='C', 0, ifelse(treatment=='P', 0, ifelse(treatment=='T', 0, ifelse(treatment=='TP', 0, 2)))),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='C', 0, ifelse(treatment=='N', 0, ifelse(treatment=='T', 0, ifelse(treatment=='TN', 0, 50)))),
         temp=ifelse(treatment=='C', 0, ifelse(treatment=='N', 0, ifelse(treatment=='P', 0, ifelse(treatment=='PN', 0, 1)))),
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()

#should N and water treatments continue all years, with the code for cessation as well so we know what they had been doing?
#cattle excluded starting in 1969 (2 years before experiment started)
esa<-read.delim("SGS_ESA.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=1, other_manipulation=0,
         n=ifelse(treatment=='control', 0, ifelse(treatment=='water', 0, 10)),
         p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='control', 0, ifelse(treatment=='N', 0, 184)),
         temp=0, precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=ifelse(calendar_year>=1976, 1, 0))%>%
  unique()

#planted in 2008 (seeds and plugs) and everything was hayed 3x per year
uk<-read.delim("SKY_UK.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=0, light=0, carbon=0, water=1, other_manipulation=1,
         n=0, p=0, k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0,
         precip=ifelse(treatment=='C', 0, ifelse(treatment=='H', 0, 30)),
         temp=ifelse(treatment=='C', 0, ifelse(treatment=='P', 0, 3)),
         precip_vari=0, precip_season=0, mow_clip=1, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=1, cessation=0)%>%
  unique()

gane<-read.delim("SVA_GANE.txt")%>%
  select(site_code, project_name, calendar_year, treatment)%>%
  mutate(community_type=0, nutrients=1, light=0, carbon=0, water=0, other_manipulation=0,
         n=ifelse(treatment=='C', 0, ifelse(treatment=='P', 0, ifelse(treatment=='LN', 0.5, ifelse(treatment=='LNP', 0.5, 5)))),
         p=ifelse(treatment=='P', 1, ifelse(treatment=='LNP', 1, ifelse(treatment=='HNP', 1, 0))),
         k=0, other_nut=0, lime=0, soil_carbon=0, CO2=0, precip=0, temp=0,
         precip_vari=0, precip_season=0, mow_clip=0, burn=0, grazed=0, fungicide=0, herb_removal=0,
         ground_light=0, other_trt=0, initial_disturbance=0, plant_mani=0, cessation=0)%>%
  unique()



###merge all datasets
combine<-rbind(bffert, bgp, biocon, bowman, clip, clonal, culardoch, e001, e002, e6, esa, events, exp1, face, fireplots,
               gane, gap2, gb, gce, herbwood, imagine, interaction, irg, kgfert, lind, mat2, megarich, mnt, mwatfer, nfert, nin, nme,
               nsfc, oface, pennings, pplots, pq, ramps, rhps, rmapc, snfert, snow, study119, study278, t7, tide, uk, wapaclip,
               warmnut, watering, wenndex, wet, yu)

write.csv(combine, "C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\sp text filesExperimentInfo_11202015.csv")









