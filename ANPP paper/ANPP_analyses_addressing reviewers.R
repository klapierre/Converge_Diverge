library(tidyverse)
library(gridExtra)
library(codyn)
library(rsq)
library(gtable)
library(grid)
#library(MASS)#loading this package disables select in tidyverse.

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")
setwd("C:\\Users\\mavolio2\\Dropbox\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw(12)) 

#read in data

# get data ----------------------------------------------------------------

#all ANPP data
anpp<-read.csv("ANPP_Oct2017_2.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

#binning the treatments into categories
trtint<-read.csv('treatment interactions_ANPP_datasets_using.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type7, trt_type5, trt_type6, trt_type)
#linking plots to treatments
anpp_expInfo<-read.csv("ExperimentInformation_ANPP_dec2017.csv")%>%
  select(-X)

Nlevels<-anpp_expInfo%>%
  select(site_code, project_name, community_type, n)%>%
  unique()%>%
  filter(n!=0)

#getting site attributes MAP, MAT and rrich
site_info<-read.csv("SiteExperimentDetails_March2019.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, site_code, project_name, community_type, MAP, MAT, rrich)

#yearly precip data
precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)

precip<-read.csv("C:\\Users\\mavolio2\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)

precip<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)


# Clean up anpp data and make calculations --------------------------------------------------------


###select the data to use

#for CDR e001/e002 selecting treatments , 6, 8, 9. For BGP dropping mowing treatments 
##all outliers were checked and are actual data.

dat2<-merge(anpp_expInfo, anpp, by=c("site_code","project_name","community_type","treatment"))%>%
  select(-nutrients, -light, -carbon, -water, -other_manipulation, -max_trt, -public, -factorial, -block)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  mutate(delete=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7|site_code=="CDR"&treatment=="2_f_u_n"|site_code=="CDR"&treatment=="3_f_u_n"|site_code=="CDR"&treatment=="4_f_u_n"|site_code=="CDR"&treatment=="5_f_u_n"|site_code=="CDR"&treatment=="7_f_u_n"|project_name=="BGP"&treatment=="u_m_c"|project_name=="BGP"&treatment=="u_m_b"|project_name=="BGP"&treatment=="u_m_n"|project_name=="BGP"&treatment=="u_m_p"|project_name=="BGP"&treatment=="b_m_c"|project_name=="BGP"&treatment=="b_m_b"|project_name=="BGP"&treatment=="b_m_n"|project_name=="BGP"&treatment=="b_m_p"|project_name=="RHPs"&calendar_year==2003, 1, 0))%>%
  filter(delete!=1)

##NOTE KBS tilling treatments did not start until 1990, 2 years after the start of the N additions and control data.
# kbs<-dat2%>%
#   filter(site_code=="KBS")%>%
#   group_by(treatment, calendar_year)%>%
#   summarise(anpp=mean(anpp))%>%
#   ungroup%>%
#   group_by(treatment)%>%
#   summarize(n=length(calendar_year))

###fix sev data
nosev<-dat2%>%
  filter(site_project_comm!="SEV_Nfert_0")

sev<-dat2%>%
  filter(site_project_comm=="SEV_Nfert_0")%>%
  select(-treatment_year)%>%
  mutate(treatment_year= ifelse(calendar_year==2004, 10, ifelse(calendar_year==2005, 11, ifelse(calendar_year==2006, 12, ifelse(calendar_year==2007, 13, ifelse(calendar_year==2008, 14, ifelse(calendar_year==2009, 15, ifelse(calendar_year==2010, 16, ifelse(calendar_year==2011, 17, 18)))))))))

all_anpp_dat<-rbind(sev, nosev)

#Calculate temporal CV
anpp_temp_cv<-all_anpp_dat%>%
  ungroup()%>%
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id, n, precip, temp, p, k)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_cv=(anpp_temp_sd/anpp_temp_mean)*100)

##calculating PD ANPP for each year
meandat<-all_anpp_dat%>%
  group_by(site_project_comm, plot_mani, treatment, treatment_year, calendar_year)%>%
  summarize(manpp=mean(anpp))

mcontrol<-meandat%>%
  filter(plot_mani==0)%>%
  mutate(contanpp=manpp)%>%
  select(-manpp)

mtrt<-meandat%>%
  filter(plot_mani!=0)

PD_anpp_yr<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(PD=((manpp-contanpp)/contanpp)*100,
         Diff=manpp-contanpp)%>%
  mutate(treatment=treatment.x)%>%
  left_join(site_info)

###Calculating PD of ANPP and CV of ANPP for each treatment plot
Means<-anpp_temp_cv%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani, n, precip, temp, p, k)%>%
  summarize(anpp_temp_cv=mean(anpp_temp_cv), 
            anpp_temp_mean=mean(anpp_temp_mean),
            anpp_temp_sd = mean(anpp_temp_sd))%>%
  ungroup()%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))

# trt_temp<-anpp_temp_cv%>%
#   filter(plot_mani>0)
# 
# CT_comp_plot<-cont_temp%>%
#   left_join(trt_temp)%>%
#   mutate(id=paste(site_code, project_name, community_type, sep="_"))%>%
#   left_join(trtint)%>%
#   mutate(PD_CV=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100,
#          PD_sd=((anpp_temp_sd-cont_temp_sd)/cont_temp_sd)*100,
#          PD_mean=((anpp_temp_mean-cont_temp_mean)/cont_temp_mean)*100)
# 

cont_temp<-Means%>%
  filter(plot_mani==0)%>%
  rename(cont_temp_cv=anpp_temp_cv,
         cont_temp_sd=anpp_temp_sd,
         cont_temp_mean=anpp_temp_mean)%>%
  select(-plot_mani, -treatment)


CT_comp<-Means%>%
  filter(plot_mani!=0)%>%
  left_join(cont_temp)%>%
  left_join(trtint)%>%
  mutate(PD_CV=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100,
         PD_sd=((anpp_temp_sd-cont_temp_sd)/cont_temp_sd)*100,
         PD_mean=((anpp_temp_mean-cont_temp_mean)/cont_temp_mean)*100)

# Getting site characteristics --------------------------------------------


##getting average production and precip vari
ave_prod<-all_anpp_dat%>%
  filter(plot_mani==0)%>%
  group_by(site_code, project_name, community_type, calendar_year)%>%
  summarize(anpp = mean(anpp))%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(manpp = mean(anpp))

precip_vari<-merge(all_anpp_dat, precip, by=c("site_code","calendar_year"))%>%
  select(site_code, project_name, community_type, calendar_year, precip_mm)%>%
  unique()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(sdppt=sd(precip_mm))

###calculate community data
#get evenness for each plot
#To calcualte evenness
#getting evenness
anpp_spc<-all_anpp_dat%>%
  select(site_project_comm)%>%
  unique()

#read in community data
community<-read.csv("SpeciesRelativeAbundance_March2019.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  right_join(anpp_spc)
spc<-unique(community$site_project_comm)
rich_even<-data.frame()

for (i in 1:length(spc)){
  subset<-community%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  rich_even<-rbind(rich_even, out)
}

trt<-community%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()
plot_mani<-all_anpp_dat%>%
  select(site_project_comm, treatment, plot_mani)%>%
  unique()

ave_even<-rich_even%>%
  left_join(trt)%>%
  group_by(site_project_comm, plot_id, treatment)%>%
  summarize(Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(site_project_comm, treatment)%>%
  summarize(Evar = mean(Evar, na.rm=T))%>%
  right_join(plot_mani)

cont_even<-ave_even%>%
  filter(plot_mani==0)%>%
  select(-plot_mani, -treatment)

site_char<-site_info%>%
  right_join(cont_even)%>%
  left_join(ave_prod)%>%
  left_join(precip_vari)

# correlating trt amount with reponses ------------------------------------
CT_comp2<-CT_comp%>%
  group_by(site_project_comm, treatment)%>%
  mutate(multnuts=sum(n, p, k))

#nitrogen - not sig for mean
with(subset(CT_comp2, trt_type6=="Nitrogen"), cor.test(n, PD_mean))
with(subset(CT_comp2, trt_type6=="Nitrogen"), plot(n, PD_mean))
#not sig for CV
with(subset(CT_comp2, trt_type6=="Nitrogen"), cor.test(n, PD_CV))
with(subset(CT_comp2, trt_type6=="Nitrogen"), plot(n, PD_CV))

#precip - not sig for mean
with(subset(CT_comp2, trt_type6=="Water"), cor.test(precip, PD_mean))
with(subset(CT_comp2, trt_type6=="Water"), plot(precip, PD_mean))
#not sig for CV
with(subset(CT_comp2, trt_type6=="Water"), cor.test(precip, PD_CV))
with(subset(CT_comp2, trt_type6=="Water"), plot(precip, PD_CV))

#temp - all only 1 C increase, but lots of vari
with(subset(CT_comp2, trt_type6=="Heat"), cor.test(temp, PD_mean))
with(subset(CT_comp2, trt_type6=="Heat"), plot(temp, PD_mean))

#mult nuts - very sig for mean
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), cor.test(multnuts, PD_mean))
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), plot(multnuts, PD_mean))
## not sig for CV
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), cor.test(multnuts, PD_CV))
with(subset(CT_comp2, trt_type6=="Multiple Nutrients"), plot(multnuts, PD_CV))
