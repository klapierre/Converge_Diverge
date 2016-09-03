#Kim
setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\clean datasets - please do not touch\\sp text files")

#Meghan
setwd("~/Dropbox/converge_diverge/datasets/FINAL_SEPT2014/clean datasets - please do not touch/sp text files")
dir()

library(gtools)
library(reshape2)
library(tidyr)
library(dplyr)

# nov 20, 2015 -checked all plots have recorded species, so the fitler abundance !=0 step will not remove any plots.

watering<-read.delim("ANG_watering.txt")%>%
  gather(species_code, abundance, sp1:sp43)%>%
  mutate(community_type=0,
         block=0)%>%
  select(-data_type)
watering_names<-read.delim("ANG_watering_specieslist.txt")
watering2<-merge(watering,watering_names, by="species_code", all=T)%>%
  select(-species_code)%>%
  filter(abundance!=0)#this drops 3 plotsin 2011 and 2013 (plots 3, 4, and 7) which had no pin hits but did have species

mat2<-read.delim("ARC_mat2.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
mat2_names<-read.delim("ARC_mat2_specieslist.txt")
mat22<-merge(mat2, mat2_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

mnt<-read.delim("ARC_mnt.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -plot_mani, -data_type, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
mnt_names<-read.delim("ARC_mnt_specieslist.txt")
mnt2<-merge(mnt, mnt_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

clonal<-read.delim("ASGA_Clonal.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -plant_mani, -plot_id1, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
clonal_names<-read.delim("ASGA_Clonal_specieslist.txt")
clonal2<-merge(clonal, clonal_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)%>%
  filter(genus_species!="No plants present",
         genus_species!="unknown grass",
         genus_species!="unknown dicot",
         genus_species!="Standing Dead biomass from current year production",
         genus_species!="unknown Asteraceae",
         genus_species!="unknonw Polygonaceae")

exp1<-read.delim("ASGA_Exp1.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -burn, -clip, -precip, -p, -dist, -patchiness, -plant_mani, -plot_id1, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp220)%>%
  mutate(community_type=0)
exp1_names<-read.delim("ASGA_Exp1_specieslist.txt")
exp12<-merge(exp1, exp1_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

nitphos<-read.csv("AZI_NitPhos.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)
  
lind<-read.delim("BAY_LIND.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -precip,-plant_mani, -plot_id1, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp56)%>%
  mutate(community_type=0)
lind_names<-read.delim("BAY_LIND_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
lind2<-merge(lind, lind_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

events<-read.delim("Bt_EVENT2.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -clip,-precip_vari, -precip_vari_season, -true_plot_mani, -plot_id1, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp56)%>%
  mutate(community_type=0)
events_names<-read.delim("Bt_EVENT2_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
events2<-merge(events, events_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

pq<-read.delim("BUX_PQ.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -true_num_manipulations, -experiment_year, -clip, -precip, -temp, -true_plot_mani, -plot_id1, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp66)%>%
  mutate(community_type=0)
pq_names<-read.delim("BUX_PQ_specieslist.txt")
pq2<-merge(pq, pq_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

pennings<-read.delim("CAR_Pennings.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -k,  -plot_id1, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(block=0)
pennings_names<-read.delim("CAR_Pennings_specieslist.txt")
pennings2<-merge(pennings, pennings_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

rmapc<-read.delim("CAU_RMAPC.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -lime,  -plot_id1, -precip, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(block=0)
rmapc_names<-read.delim("CAU_RMAPC_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
rmapc2<-merge(rmapc, rmapc_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

biocon<-read.delim("CDR_biocon.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -c, -plant_mani, -true_plot_mani, -plot_mani, -data_type, -species_num, -true_num_manipulations)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(block=0,
         community_type=0)
biocon_names<-read.delim("CDR_biocon_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
biocon2<-merge(biocon, biocon_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

e001<-read.csv("CDR_e001.csv")%>%
  mutate(block=0)%>%
  select(-data_type)%>%
  filter(abundance!=0)
e001_names<-read.csv("CDR_e001_e002_specieslist.csv")
e0012<-merge(e001, e001_names, by="spcode")%>%
  filter(abundance!=0)%>%
  select(-spcode)%>%
  filter(genus_species!="Forb seedlings",
         genus_species!="Fungi ",
         genus_species!="Miscellaneous grasses",
         genus_species!="Miscellaneous herbs",
         genus_species!="Miscellaneous legumes",
         genus_species!="Miscellaneous litter",
         genus_species!="Miscellaneous rushes",
         genus_species!="Miscellaneous sedges",
         genus_species!="Miscellaneous sp.",
         genus_species!="Miscellaneous woody")

e002<-read.delim("CDR_e002.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -true_num_manipulations, -experiment_year, -p, -k, -lime, -n, -other_nut, -burn, -herb_removal, -true_plot_mani, -plot_mani, -cessation, -dist, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(block=0)
e002_names<-read.delim("CDR_e002_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
e0022<-merge(e002, e002_names,by="species_code", all=T)%>%
  filter(abundance!=0, calendar_year<1992)%>%##drops everything once cessation starts
  mutate(spcode=genus_species)%>%
  select(-species_code, -genus_species)
e0023<-merge(e0022, e001_names, by="spcode")%>%
  select(-spcode)

megarich<-read.delim("CEH_Megarich.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -true_num_manipulations, -experiment_year, -clip, -c, -temp, -n, -p, -k, -true_plot_mani, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(block=0,
         community_type=0)
megarich_names<-read.delim("CEH_Megarich_specieslist.txt")
megarich2<-merge(megarich, megarich_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

imagine<-read.delim("CLE_imagine.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -c, -data_type,-plot_mani, -plant_mani, -precip, -temp, -species_num)%>%
  gather(species_code, abundance, sp1:sp12)%>%
  mutate(community_type=0)
imagine_names<-read.delim("CLE_imagine_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
imagine2<-merge(imagine, imagine_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

culardoch<-read.delim("CUL_culardoch.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -burn, -clip, -n, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp43)%>%
  mutate(community_type=0)
culardoch_names<-read.delim("CUL_culardoch_specieslist.txt")
culardoch2<-merge(culardoch, culardoch_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

gap2<-read.delim("DCGS_gap2.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -l, -data_type,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp125)%>%
  mutate(community_type=0)
gap2_names<-read.delim("DCGS_gap2_specieslist.txt")
gap22<-merge(gap2, gap2_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

nsfc<-read.delim("DL_NSFC.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -precip, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp53)%>%
  mutate(community_type=0)
nsfc_names<-read.delim("DL_NSFC_specieslist.txt")%>%
  mutate(species_code=tolower(species_list))
nsfc2<-merge(nsfc, nsfc_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code, -species_list)

warmnut<-read.delim("Finse_WarmNut.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -k, -temp, -data_type,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp228)%>%
  mutate(community_type=0)
warmnut_names<-read.delim("Finse_WarmNut_specieslist.txt")
warmnut2<-merge(warmnut, warmnut_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

face<-read.delim("GVN_FACE.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -c, -data_type,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0,
         block=0)
face_names<-read.delim("GVN_FACE_specieslist.txt")
face2<-merge(face, face_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

nde<-read.csv("IMGERS_NDE.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0)%>%
  filter(abundance!=0)

yu<-read.delim("IMGERS_Yu.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -k, -plot_mani, -data_type, -species_num)%>%
  gather(species_code, abundance, sp1:sp29)%>%
  mutate(community_type=0)
yu_names<-read.delim("IMGERS_Yu_specieslist.txt")
yu2<-merge(yu, yu_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

study119<-read.delim("JRN_Study119.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -data_type,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)
study119_names<-read.delim("JRN_Study119_specieslist.txt")
study1192<-merge(study119, study119_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

study278<-read.delim("JRN_study278.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -precip, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
study278_names<-read.delim("JRN_study278_specieslist.txt")
study2782<-merge(study278, study278_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

gce<-read.delim("JSP_GCE2.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -precip, -temp, -c, -data_type,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp92)%>%
  mutate(community_type=0, block=0)
gce_names<-read.delim("JSP_GCE2_specieslist.txt")
gce2<-merge(gce, gce_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

wapaclip<-read.delim("KAEFS_WaPaClip.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -clip, -precip, -temp, -data_type,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
wapaclip_names<-read.delim("KAEFS_WaPaClip_specieslist.txt")
wapaclip2<-merge(wapaclip, wapaclip_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

t7<-read.delim("KBS_T7.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -dist,-n, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
t7_names<-read.delim("KBS_T7_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
t72<-merge(t7, t7_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

bffert<-read.delim("KLU_BFFert.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -p,-n,-k,-herb_removal,-plot_id1, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
bffert_names<-read.delim("KLU_BFFert_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
bffert2<-merge(bffert, bffert_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

kgfert<-read.delim("KLU_KGFert.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -p,-n,-k,-fungicide, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
kgfert_names<-read.delim("KLU_KGFert_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
kgfert2<-merge(kgfert, kgfert_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

bgp<-read.delim("KNZ_BGP.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -p,-n,-burn,-clip, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp220)%>%
  mutate(community_type=0,
         block=0)
bgp_names<-read.delim("KNZ_BGP_specieslist.txt")
bgp2<-merge(bgp, bgp_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

irg<-read.delim("KNZ_IRG.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -precip, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp220)%>%
  mutate(block=0)
irg_names<-read.delim("KNZ_IRG_specieslist.txt")
irg2<-merge(irg, irg_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

gfp<-read.csv("KNZ_KNP_GFP.csv")%>%
  select(-data_type)%>%
  mutate(block=0)%>%
  filter(abundance!=0)

pplots<-read.csv("KNZ_PPLOTS.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)

ramps<-read.csv("KNZ_Ramps.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)

rhps<-read.delim("KNZ_RHPs.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -soil_depth, -data_type,-plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
rhps_names<-read.delim("KNZ_RHPs_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
rhps2<-merge(rhps, rhps_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

e6<-read.delim("KUFS_E6.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p,-plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp168)
e6_names<-read.delim("KUFS_E6_specieslist.txt")
e62<-merge(e6, e6_names,by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

clip<-read.delim("LATNJA_CLIP.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -temp, -data_type, -plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp227)
clip_names<-read.delim("LATNJA_CLIP_specieslist.txt")
clip2<-merge(clip, clip_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

pme<-read.csv("LEFT_PME.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0)%>%
  filter(abundance!=0)

herbwood<-read.delim("LG_HerbWood.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n,-precip, -p, -k, -data_type, -plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(block=0,
         community_type=0)
herbwood_names<-read.delim("LG_HerbWood_specieslist.txt")
herbwood2<-merge(herbwood, herbwood_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

fireplots<-read.delim("MAERC_fireplots.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n,-burn, -p, -clip, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
fireplots_names<-read.delim("MAERC_fireplots_specieslist.txt")
fireplots2<-merge(fireplots, fireplots_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

mwatfer<-read.csv("MNR_watfer.csv")%>%
  mutate(genus_species=species_name)%>%
  select(-species_name, -species)%>%
  filter(abundance!=0)

wet<-read.delim("NANT_wet.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -data_type, -plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp101)%>%
  mutate(block=0)
wet_names<-read.delim("NANT_wet_specieslist.txt")
wet2<-merge(wet, wet_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

gb<-read.delim("NGBER_gb.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -precip_vari_season, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp21)%>%
  mutate(community_type=0)
gb_names<-read.delim("NGBER_gb_specieslist.txt")
gb2<-merge(gb, gb_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

herbdiv<-read.csv("NIN_herbdiv.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)

ccd<-read.delim("NTG_CCD.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -clip,-precip, -temp, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp94)%>%
  mutate(community_type=0,
         block=0)
ccd_names<-read.delim("NTG_CCD_specieslist.txt")
ccd2<-merge(ccd, ccd_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

nfert<-read.delim("NWT_246NFert.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n,  -data_type, -plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
nfert_names<-read.delim("NWT_246NFert_specieslist.txt")
nfert2<-merge(nfert, nfert_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

bowman<-read.delim("NWT_bowman.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -p, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)
bowman_names<-read.delim("NWT_bowman_specieslist.txt")
bowman2<-merge(bowman, bowman_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

snow<-read.delim("NWT_snow.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -temp, -precip,-data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
snow_names<-read.delim("NWT_snow_specieslist.txt")
snow2<-merge(snow, snow_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

oface<-read.delim("ORNL_FACE.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -c, -data_type, -plot_mani, -species_num, -plot_id1)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0, block=0)
oface_names<-read.delim("ORNL_FACE_specieslist.txt")
oface2<-merge(oface, oface_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

tide<-read.delim("PIE_Tide.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0)
tide_names<-read.delim("PIE_Tide_specieslist.txt")
tide2<-merge(tide, tide_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

interaction<-read.delim("RIO_interaction.txt")%>%
  select(-n, -precip, -precip_vari, -plot_mani, -data_type)%>%
  gather(species_code, abundance, sp1:sp10)%>%
  mutate(block=0)
interaction_names<-read.delim("RIO_interaction_specieslist.txt")
interaction2<-merge(interaction, interaction_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)  

lucero<-read.csv("SCL_Lucero.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0)%>%
  filter(abundance!=0)%>%
  filter(genus_species!="Standing.dead")

ter<-read.csv("SCL_TER.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0)%>%
  filter(abundance!=0)

cxn<-read.csv("SERC_CXN.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)

tmece<-read.csv("SERC_TMECE.csv")%>%
  select(-data_type)%>%
  mutate(block=0)%>%
  filter(abundance!=0)

snfert<-read.delim("SEV_NFert.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0,
         block=0)
snfert_names<-read.delim("SEV_NFert_specieslist.txt")%>%
  mutate(species_code=tolower(species_code))
snfert2<-merge(snfert, snfert_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)
  
wenndex<-read.delim("SEV_WENNDEx.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -experiment_year, -n, -temp, -precip, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0,
         block=0)
wenndex_names<-read.delim("SEV_WENNDEx_specieslist.txt")
wenndex2<-merge(wenndex, wenndex_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

grazeprecip<-read.csv("SFREC_GrazePrecip.csv")%>%
  select(-data_type)%>%
  filter(abundance!=0)

uk<-read.delim("SKY_UK.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulation, -num_manipulations, -temp, -precip, -plot_id1, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp26)%>%
  mutate(community_type=0)
uk_names<-read.delim("SKY_UK_specieslist.txt")
uk2<-merge(uk, uk_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

nitrogen<-read.csv("SR_Nitrogen.csv")%>%
  select(-data_type)%>%
  filter(abundance!=0, genus_species!="UNKNOWN SPECIES", genus_species!="GRASS SPECIES", genus_species!="FORB SPECIES")

water<-read.csv("SR_WATER.csv")%>%
  select(-data_type)%>%
  filter(abundance!=0, genus_species!="UNKNOWN SPECIES", genus_species!="GRASS SPECIES", genus_species!="FORB SPECIES")

gane<-read.delim("SVA_GANE.txt")%>%
  select(-id, -nutrients, -light, -carbon, -water, -other_manipulations, -num_manipulations, -experiment_year, -n, -p, -data_type, -plot_mani, -species_num)%>%
  gather(species_code, abundance, sp1:sp232)%>%
  mutate(community_type=0,
         block=0)
gane_names<-read.delim("SVA_GANE_specieslist.txt")
gane2<-merge(gane, gane_names, by="species_code", all=T)%>%
  filter(abundance!=0)%>%
  select(-species_code)

tface<-read.csv("TAS_FACE.csv")%>%
  select(-data_type, -X)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)

lovegrass<-read.csv("TRA_Lovegrass.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0)%>%
  filter(abundance!=0)

nitadd<-read.csv("YMN_NitAdd.csv")%>%
  select(-data_type)%>%
  mutate(community_type=0, block=0)%>%
  filter(abundance!=0)

#merge all datasets
combine<-rbind(bffert2, bgp2, biocon2, bowman2, ccd2, clip2, clonal2, culardoch2, cxn, e0012, e0023, e62, events2, exp12, face2, fireplots2, gane2, gap22, gb2, gce2, gfp, grazeprecip, herbdiv, herbwood2, imagine2, interaction2, irg2, kgfert2, lind2, lovegrass, lucero, mat22, megarich2, mnt2, mwatfer, nde, nfert2, nitadd, nitphos, nitrogen, nsfc2, oface2,pennings2, pme, pplots, pq2, ramps, rhps2, rmapc2, snfert2, snow2, study1192, study2782, t72, ter, tface, tide2, tmece, uk2 ,wapaclip2, warmnut2, watering2, water, wenndex2, wet2, yu2)%>%
  mutate(genus_species=tolower(genus_species))

#take2<-aggregate(abundance~site_code+project_name+community_type, sum, data=combine)

write.csv(combine, "C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\SpeciesRawAbundance_Sept2016.csv")

###get species list
species_list<-combine%>%
  select(site_code, project_name, genus_species)%>%
  unique()

write.csv(species_list, "C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\SpeciesList_Sept2016.csv")

###Getting Relative Cover
totcov<-combine%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, treatment, block, plot_id)%>%
  summarise(totcov=sum(abundance))

relcov<-merge(totcov, combine, by=c("site_code", "project_name", "community_type", "calendar_year", "treatment_year", "treatment", "block", "plot_id"), all=T)%>%
  mutate(relcov=abundance/totcov)%>%
  select(-abundance, -totcov)

write.csv(relcov, "C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Sept2016.csv")

##for Codyn dataset
expinfo<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/ExperimentInformation_Nov2015.csv")
forcodyn<-merge(expinfo, relcov, by=c("site_code","project_name","community_type","treatment_year","calendar_year","treatment"))%>%
  select(site_code,project_name,community_type,calendar_year,treatment_year,genus_species,relcov,plot_mani, plot_id)

##get sp1 sp2
species<-forcodyn%>%
  select(site_code,project_name, genus_species)%>%
  mutate(exp=paste(site_code, project_name, sep="::"))%>%
  unique()

exp<-species%>%
  select(site_code,project_name)%>%
  mutate(exp=paste(site_code, project_name, sep="::"))%>%
  unique()

forcodynsp<-data.frame(row.names = 1)

for (i in 1: length(exp$exp)) {
  
  subset<-species%>%
    filter(exp==exp$exp[i])
  
  newnames<-subset%>%
    mutate(sp="sp", num=seq(1,nrow(subset), by=1), species=paste(sp, num, sep=""))
  
  forcodynsp<-rbind(newnames, forcodynsp)
}

forcodyn2<-merge(forcodyn, forcodynsp, by=c("site_code", "project_name", "genus_species"))%>%
  select(-genus_species, -sp, -num, -exp)%>%
  #there are some repeat species!GAH
  tbl_df()%>%
  group_by(site_code, project_name, community_type, calendar_year, treatment_year, plot_mani, plot_id, species)%>%
  summarize(relcov=mean(relcov))

write.csv(forcodyn2,"~/Dropbox/CoDyn/R files/11_06_2015_v7/corre_relcov.csv")
