setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\clean datasets - please do not touch\\anpp text files")

library(gtools)
library(reshape2)
library(tidyr)
library(dplyr)

nme<-read.delim("Alps_NME_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
watering<-read.delim("ANG_watering_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
clonal<-read.delim("ASGA_Clonal_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
exp1<-read.delim("ASGA_Exp1_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
lind<-read.delim("BAY_LIND_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
events<-read.delim("Bt_EVENT2_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
rmapc<-read.delim("CAU_RMAPC_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
biocon<-read.delim("CDR_biocon_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
e001<-read.delim("CDR_e001_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
e002<-read.delim("CDR_e002_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
megarich<-read.delim("CEH_Megarich_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
imgers<-read.delim("IMGERS_Yu_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
imagine<-read.delim("CLE_imagine_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
nsfc<-read.delim("DL_NSFC_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
wapaclip<-read.delim("KAEFS_WaPaClip_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
t7<-read.delim("KBS_T7_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
kgfert<-read.delim("KLUG_KGFert_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
bgp<-read.delim("KNZ_BGP_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
irg<-read.delim("KNZ_IRG_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
pplots<-read.delim("KNZ_PPLOTS_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
ramps<-read.delim("KNZ_Ramps_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
rhps<-read.delim("KNZ_RHPs_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
nfert<-read.delim("NWT_246NFert_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
bowman<-read.delim("NWT_bowman_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
oface<-read.delim("ORNL_FACE_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
tide<-read.delim("PIE_Tide_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
snfert<-read.delim("SEV_NFert_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
wenndex<-read.delim("SEV_WENNDEx_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
uk<-read.delim("SKY_UK_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
fireplots<-read.delim("MAERC_fireplots_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
wet<-read.delim("NANT_wet_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
gb<-read.delim("NGBER_gb_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)
rio<-read.delim("RIO_interaction_anpp.txt")%>%
  select(site_code, project_name, community_type, treatment_year, calendar_year, plot_id, anpp)
mnr<-read.delim("MNR_watfer_anpp.txt")%>%
  select(site_code, project_name, treatment_year, calendar_year, plot_id, anpp)%>%
  mutate(community_type=0)

e001$community_type<-as.character(e001$community_type)
clonal$plot_id1<-as.character(clonal$plot_id1)
events$precip_vari_season<-as.character(events$precip_vari_season)
watering$precip_season<-as.character(watering$precip_season)
biocon$block<-as.character(biocon$block)

combine_anpp<-smartbind(bgp, biocon, bowman, clonal, e001, e002, events, exp1, imagine, irg, kgfert, lind, megarich, nfert, nme, nsfc, oface, pplots, qiang, ramps, rhps, rmapc, snfert, t7, tide, uk, wapaclip, watering, wenndex, fill=0)

write.csv(combine_anpp, "AllAnppData_02172015.csv")

# take2<-aggregate(light~site_code+project_name, sum, data=take1)

# str(bgp)
# str(biocon)
# str(bowman)
# str(clonal)
# str(e001)
# str(e002)
# str(events)
# str(exp1)
# str(face)
# str(imagine)
# str(irg)
# str(kgfert)
# str(lind)
# str(megarich)
# str(nfert)
# str(nsfc)
# str(oface)
# str(pplots)
# str(qiang)
# str(ramps)
# str(rhps)
# str(snfert)
# str(t7)
# str(tide)
# str(uk)
# str(wapaclip)
# str(watering)
# str(wenndex)
