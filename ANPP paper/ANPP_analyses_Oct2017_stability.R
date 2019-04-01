library(tidyverse)
library(ggplot2)
library(gridExtra)
library(gtable)
library(devtools)
install_github("NCEAS/codyn", ref = github_pull(83))
library(codyn)
library(lme4)
# library(gtools)
# library(grid)
library(lmerTest)
# library(MuMIn)
# library(ppcor)

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw(14)) 

#read in data

# get data ----------------------------------------------------------------


anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)

site_info<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, MAP, MAT, rrich)

anpp<-read.csv("ANPP_Oct2017_2.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

trtint<-read.csv('treatment interactions_ANPP_datasets_using.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type7, trt_type5, trt_type6, trt_type)

#no longer using prism data
# precip<-read.csv('~/Dropbox/converge_diverge/datasets/LongForm/climate/ANPP_PrecipData.csv')
# 
# precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\ANPP_PrecipData.csv")%>%
#   mutate(site_code=ï..site_code)%>%
#   select(-ï..site_code)

precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)

precip<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)

# clean up anp data --------------------------------------------------------


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

nosev<-dat2%>%
  filter(site_project_comm!="SEV_Nfert_0")

sev<-dat2%>%
  filter(site_project_comm=="SEV_Nfert_0")%>%
  select(-treatment_year)%>%
  mutate(treatment_year= ifelse(calendar_year==2004, 10, ifelse(calendar_year==2005, 11, ifelse(calendar_year==2006, 12, ifelse(calendar_year==2007, 13, ifelse(calendar_year==2008, 14, ifelse(calendar_year==2009, 15, ifelse(calendar_year==2010, 16, ifelse(calendar_year==2011, 17, 18)))))))))

all_anpp_dat<-rbind(sev, nosev)

sites<-all_anpp_dat%>%
  group_by(site_project_comm,site_code, calendar_year, treatment, plot_mani)%>%
  summarize(anpp=mean(anpp))%>%
  mutate(spc_trt=paste(site_project_comm, treatment, sep="::"))%>%
  group_by(site_project_comm, site_code, treatment)%>%
  summarize(len=length(calendar_year))
 

# ggplot(data=all_anpp_dat, aes(anpp))+
#   geom_histogram()+
#   facet_wrap(~site_project_comm, ncol=4, scales="free")


#write.csv(all_anpp_dat, "ANPP_6yrs_Dec2017.csv")
# 
# anpp_trts<-all_anpp_dat%>%
#   select(site_project_comm, treatment)%>%
#   unique

# calculate temporal spatial cv and effect sizes ---------------------------------------

#calculate spatail
anpp_spatial<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_mani)%>%
  summarize(anpp_sp_sd=sd(anpp, na.rm=T),
            anpp_sp_mean=mean(anpp, na.rm=T),
            anpp_sp_stab=anpp_sp_mean/anpp_sp_sd)%>%
  select(-anpp_sp_sd, -anpp_sp_mean)%>%
  filter(treatment_year!=0)

lastyr<-anpp_spatial%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  group_by(site_project_comm)%>%
  mutate(maxyr=max(calendar_year))%>%
  filter(maxyr==calendar_year)

other_year<-anpp_spatial%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  group_by(site_project_comm)%>%
  mutate(otheryr=max(calendar_year)-4)%>%
  filter(otheryr==calendar_year)


#Calculate temporal
anpp_temp_cv<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_stab=anpp_temp_mean/anpp_temp_sd)%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarize(anpp_temp_stab=mean(anpp_temp_stab, na.rm=T))

##calculating effect sizes
meandat<-all_anpp_dat%>%
  group_by(site_project_comm, plot_mani, treatment, treatment_year, calendar_year)%>%
  summarize(manpp=mean(anpp))

mcontrol<-meandat%>%
  filter(plot_mani==0)%>%
  mutate(contanpp=manpp)%>%
  select(-manpp)

mtrt<-meandat%>%
  filter(plot_mani!=0)

logRR<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(logrr=abs(log(manpp/contanpp)),
         RR=((manpp-contanpp)/contanpp)*100)%>%
  group_by(site_project_comm, treatment.x)%>%
  summarise(mlogrr=mean(logrr),
            mrr=sd(RR))%>%
  mutate(treatment=treatment.x)%>%
  select(-treatment.x)

logRRsp<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(logrr=abs(log(manpp/contanpp)))%>%
  group_by(site_project_comm, treatment.x, calendar_year)%>%
  summarise(mlogrr=mean(logrr),
            sdlogrr=sd(logrr))%>%
  mutate(treatment=treatment.x)%>%
  ungroup()%>%
  select(-treatment.x)


# overall effect of vari --------------------------------------------------
cont_temp<-anpp_temp_cv%>%
  filter(plot_mani==0)%>%
  mutate(cont_temp_stab=anpp_temp_stab)%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, cont_temp_stab)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(-treatment)

trt_temp<-anpp_temp_cv%>%
  filter(plot_mani>0)

tograph1_temp<-merge(cont_temp, trt_temp, by=c("site_code", 'project_name',"community_type"))%>%
  mutate(id=paste(site_code, project_name, community_type, sep="_"))

tograph_temp<-merge(tograph1_temp, trtint, by=c("site_project_comm","treatment"))

##spatial

cont_spat<-lastyr%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_stab=anpp_sp_stab)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_stab, calendar_year)%>%
  select(-treatment)

trt_spat<-lastyr%>%
  filter(plot_mani>0)

tograph1_spat<-merge(cont_spat, trt_spat, by=c("site_project_comm","calendar_year"))

tograph_spat<-merge(tograph1_spat, trtint, by=c("site_project_comm","treatment"))
###other year for spatial
cont_spatother<-other_year%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_cv=anpp_sp_cv)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_cv, calendar_year)%>%
  select(-treatment)

trt_spatother<-other_year%>%
  filter(plot_mani>0)
tograph1_otherspat<-merge(cont_spatother, trt_spatother, by=c("site_project_comm","calendar_year"))

tograph_ottherspat<-merge(tograph1_otherspat, trtint, by=c("site_project_comm","treatment"))

# ####just overall what are the effects of the treatments on spatial and temporal heterogeneity?
# ##testing for differences
# spat_bar<-tograph_spat%>%
#   mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)
#
# t.test(abs(spat_bar$PC), mu=0)
# t.test(spat_bar$PC, mu=0)
#
# temp_bar<-tograph_temp%>%
#   mutate(PC=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100)
#
# t.test(abs(temp_bar$PC), mu=0)
# t.test(temp_bar$PC, mu=0)
#
# ##other year
# other_bar<-tograph_ottherspat%>%
#   mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)
#
#
# t.test(abs(other_bar$PC), mu=0)
#
# #do t-test do the PC differ from zero?
# #spatial
# irr<-subset(spat_bar, trt_type6=="Water")
# t.test(irr$PC, mu=0)
# nit<-subset(spat_bar, trt_type6=="Nitrogen")
# t.test(nit$PC, mu=0)
# nuts<-subset(spat_bar, trt_type6=="Multiple Nutrients")
# t.test(nuts$PC, mu=0)
#
#
# #temporal model
# irr<-subset(temp_bar, trt_type6=="Water")
# t.test(irr$PC, mu=0)
# nit<-subset(temp_bar, trt_type6=="Nitrogen")
# t.test(nit$PC, mu=0)
# nuts<-subset(temp_bar, trt_type6=="Multiple Nutrients")
# t.test(nuts$PC, mu=0)
#
#
# tograph_spat_bar<-tograph_spat%>%
#   mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)%>%
#   group_by(trt_type6)%>%
#   summarize(P.C=mean(PC),
#             sdd=sd(PC),
#             num=length(PC))%>%
#   mutate(se=sdd/sqrt(num))%>%
#   filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")
#
# tograph_temp_bar<-tograph_temp%>%
#   mutate(PC=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100)%>%
#   group_by(trt_type6)%>%
#   summarize(P.C=mean(PC),
#             sdd=sd(PC),
#             num=length(PC))%>%
#   mutate(se=sdd/sqrt(num))%>%
#   filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")
#
#
# #graphing this
# temp_pc<-
#   ggplot(data=tograph_temp_bar, aes(x=trt_type6, y=P.C, fill=trt_type6))+
#   geom_bar(position=position_dodge(), stat="identity")+
#   geom_errorbar(aes(ymin=P.C-se, ymax=P.C+se),position= position_dodge(0.9), width=0.2)+
#   ylab("Percent Change of Temporal Variability")+
#   theme(axis.text.x=element_text(angle=45, hjust=1))+
#   scale_fill_manual(values=c("darkred","orange","dodgerblue"))+
#   xlab("Treatment")+
#   ggtitle("Temporal")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(legend.position = "none")+
#   scale_y_continuous(limits=c(-35,100))
# spat_pc<-
#   ggplot(data=tograph_spat_bar, aes(x=trt_type6, y=P.C, fill=trt_type6))+
#   geom_bar(position=position_dodge(), stat="identity")+
#   geom_errorbar(aes(ymin=P.C-se, ymax=P.C+se),position= position_dodge(0.9), width=0.2)+
#   ylab("Percent Change of Spatial Variability")+
#   theme(axis.text.x=element_text(angle=45, hjust=1))+
#   scale_fill_manual(values=c("darkred","orange","dodgerblue"))+
#   xlab("Treatment")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   ggtitle("Spatial")+
#   theme(legend.position = "none")+
#   geom_text(x=1, y=75, label="*", size=8)+
#   geom_text(x=2, y=90, label="*", size=8)+
#   scale_y_continuous(limits=c(-35,100))
#
# grid.arrange(temp_pc, spat_pc, ncol=2)



# Q1 how does gcds effect temporal or spatial vari? -----------------------

###temporal analysis


##t-test - do the slopes differ from 1?
#Temporal control versus treatment
# do the lines differ from a slope of 1?
temp.lm<-lm(anpp_temp_stab~cont_temp_stab, data=tograph_temp)
my.slope <- summary(temp.lm)$coef["cont_temp_stab", c("Estimate", "Std. Error")]
my.df <- summary(temp.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001

#spatail control versus treatment
# do the lines differ from a slope of 1?
spat.lm<-lm(anpp_sp_stab~cont_sp_stab, data=tograph_spat)
my.slope <- summary(spat.lm)$coef["cont_sp_stab", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001

###test another year for spatial data

spat.lm<-lm(anpp_sp_stab~cont_sp_stab, data=tograph_ottherspat)
my.slope <- summary(spat.lm)$coef["cont_sp_cv", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001


###graphing this
temp<-
ggplot(data=tograph_temp, aes(x=cont_temp_stab, y=anpp_temp_stab, color = trt_type7))+
  scale_color_manual(name = "GCD Treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  geom_point(size=2)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="black", size = 2)+
  geom_smooth(data=subset(tograph_temp, trt_type6 =="Nitrogen"), method="lm", se=F, color="green3", size = 1)+
  geom_smooth(data=subset(tograph_temp, trt_type6 =="Multiple Nutrients"), method="lm", se=F, color="orange", size = 1)+
  geom_smooth(data=subset(tograph_temp, trt_type6 =="Water"), method="lm", se=F, color="blue", size = 1)+
  ylab("Temporal CV Treatment Plots")+
  xlab("Temporal CV Control Plots")+
  ggtitle("Temporal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", x = 25, y = 95, label= "A", size = 8)

spat<-
ggplot(data=tograph_spat, aes(x=cont_sp_stab, y=anpp_sp_stab, color = trt_type7))+
  scale_color_manual(name = "GCD Treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  geom_point(size=2)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="black", size = 2)+
  geom_smooth(data=subset(tograph_spat, trt_type6 =="Nitrogen"), method="lm", se=F, color="green3", size = 1)+
  geom_smooth(data=subset(tograph_spat, trt_type6 =="Multiple Nutrients"), method="lm", se=F, color="orange", size = 1)+
  geom_smooth(data=subset(tograph_spat, trt_type6 =="Water"), method="lm", se=F, color="blue", size = 1)+
  ylab("Spatial CV Treatment Plots")+
  xlab("Spatial CV Control Plots")+
  ggtitle("Spatial")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", x = 10, y = 135, label= "B", size = 8)

grid.arrange(temp, spat, ncol=1)


###looking at three well replicated treatmetns.
#nitrogen - temporal
temp.lm<-lm(anpp_temp_stab~cont_temp_stab, data=subset(tograph_temp, trt_type6=="Nitrogen"))
my.slope <- summary(temp.lm)$coef["cont_temp_stab", c("Estimate", "Std. Error")]
my.df <- summary(temp.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# no p = 0.068

# nuts temporal
temp.lm<-lm(anpp_temp_stab~cont_temp_stab, data=subset(tograph_temp, trt_type6=="Multiple Nutrients"))
my.slope <- summary(temp.lm)$coef["cont_temp_stab", c("Estimate", "Std. Error")]
my.df <- summary(temp.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p = 0.003

#water temporal
temp.lm<-lm(anpp_temp_stab~cont_temp_stab, data=subset(tograph_temp, trt_type6=="Water"))
my.slope <- summary(temp.lm)$coef["cont_temp_stab", c("Estimate", "Std. Error")]
my.df <- summary(temp.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p = 1


#nitrogen spatial
spat.lm<-lm(anpp_sp_cv~cont_sp_cv, data=subset(tograph_spat, trt_type6=="Nitrogen"))
my.slope <- summary(spat.lm)$coef["cont_sp_cv", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# no p = 0.467

#nuts
spat.lm<-lm(anpp_sp_cv~cont_sp_cv, data=subset(tograph_spat, trt_type6=="Multiple Nutrients"))
my.slope <- summary(spat.lm)$coef["cont_sp_cv", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001

#Water
spat.lm<-lm(anpp_sp_cv~cont_sp_cv, data=subset(tograph_spat, trt_type6=="Water"))
my.slope <- summary(spat.lm)$coef["cont_sp_cv", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# no p = 0.195

# tograph_temp_trt<-tograph_temp%>%
#   filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")
# 
# tograph_spat_trt<-tograph_spat%>%
#   filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")
# temp_trt<-
# ggplot(data=tograph_temp_trt, aes(x=cont_temp_cv, y=anpp_temp_cv))+
#   geom_point(size=2)+
#   geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
#   geom_smooth(method="lm", se=F, color="black")+
#   ylab("Temporal CV Treatment Plots")+
#   xlab("Temporal CV Control Plots")+
#   ggtitle("A) Temporal")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_wrap(~trt_type6)
# spat_trt<-
# ggplot(data=tograph_spat_trt, aes(x=cont_sp_cv, y=anpp_sp_cv))+
#   geom_point(size=2)+
#   geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
#   geom_smooth(method="lm", se=F, color="black")+
#   ylab("Spatial CV Treatment Plots")+
#   xlab("Spatial CV Control Plots")+
#   ggtitle("B) Spatial")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_wrap(~trt_type6)
# 
# grid.arrange(temp_trt, spat_trt, ncol=2)

###spatial through time.
# cont_spat_all<-anpp_spatial%>%
#   mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
#   filter(plot_mani==0)%>%
#   mutate(cont_sp_cv=anpp_sp_cv)%>%
#   ungroup()%>%
#   select(site_project_comm, treatment, cont_sp_cv, calendar_year, treatment_year)%>%
#   select(-treatment)
# 
# trt_spat_all<-anpp_spatial%>%
#   mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
#   filter(plot_mani>0)
# 
# tograph1_spat_all<-merge(cont_spat_all, trt_spat_all, by=c("site_project_comm","treatment_year","calendar_year"))
# 
# tograph_spat_all<-merge(tograph1_spat_all, trtint, by=c("site_project_comm","treatment"))%>%
#   mutate(logrr=log(anpp_sp_cv/cont_sp_cv))%>%
#   mutate(spc_t=paste(site_project_comm, treatment, sep=""))
# 
# ggplot(data=tograph_spat_all, aes(x=treatment_year, y=logrr))+
#   geom_line(aes(group=spc_t), size=0.5)+
#   ylab("Log RR of Spatial CV")+
#   geom_abline(slope=0, intercept=0, size=1)+
#   xlab("Treatment Year")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Q2 what is the relationship between control CV and effect size? ---------

##temporal
tograph_log1_temp<-merge(logRR, cont_temp, by="site_project_comm")
tograph_log2_temp<-merge(tograph_log1_temp, trtint, by=c("site_project_comm","treatment"))
tograph_log_temp<-merge(tograph_log2_temp, site_info, by="site_project_comm")

#test the relationship between control_temp and effect size

temp_effect <- lm(mlogrr ~ cont_temp_stab, data = tograph_log_temp)
summary(temp_effect)

# map_effect <- lm(cont_temp_cv ~ MAP, data = tograph_log_temp)
# summary(map_effect)
# 
# map_temp_effect <- lm(mlogrr ~ MAP, data = tograph_log_temp)
# summary(map_temp_effect)

# ##do partial correlation to see the correlation between cv control and logrr given MAP
# pcordata<-tograph_log_temp[,c(3,8,12)]
# with(tograph_log_temp, pcor.test(cont_temp_cv, mlogrr, MAP))
# pcor(pcordata)

#looking at three seperate GCDs

summary(lm(mlogrr ~ cont_temp_cv, 
                  data = subset(tograph_log_temp, trt_type6=="Nitrogen")))
summary(lm(mlogrr ~ cont_temp_cv, 
           data = subset(tograph_log_temp, trt_type6=="Water")))

summary(lm(mlogrr ~ cont_temp_cv, 
           data = subset(tograph_log_temp, trt_type6=="Multiple Nutrients")))


##graphing this

  ggplot(data=tograph_log_temp, aes(x=cont_temp_stab, y=mlogrr, color = trt_type7))+
  geom_point(size=2)+
  scale_color_manual(name = "GCD Treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  ylab("Log RR")+
  xlab("Temporal CV Control Plots")+
  geom_smooth(method="lm", color="black", se=F, size = 2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

tograph_log_temp_trt<-tograph_log_temp%>%
  filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")

ggplot(data=tograph_log_temp_trt, aes(x=cont_temp_cv, y=mlogrr))+
  geom_point(size=2)+
  ylab("Log RR")+
  xlab("Temporal CV Control Plots")+
  #geom_smooth(method="lm", color="black", se=F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt_type6)

# precipitation analysis --------------------------------------------------

#precip analysis #1986 in CDR has no precip data, this one year is being dropped.

#drop irrigation treatments becuase it is confusing to add water and then test against precip OR maybe just KNZ because not the same amount of water each year.

#look at average change in sensitivity by treatments.

anpp_precip<-merge(all_anpp_dat, precip, by=c("site_code","calendar_year"))%>%
  mutate(trt=ifelse(plot_mani==0,"C","T"))%>%
  group_by(site_project_comm, trt, calendar_year, precip_mm, treatment, plot_mani)%>%
  summarize(anpp=mean(anpp))%>%
  mutate(spc_trt=paste(site_project_comm, treatment, sep="::"))

ggplot(data=anpp_precip, aes(x=precip_mm, y=anpp, group=treatment, color=trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+
  facet_wrap(~site_project_comm, ncol=8, scales = "free")


##get slopes for each treatment including controls
spc<-unique(anpp_precip$spc_trt)
lm.slopes<-data.frame()
for (i in 1:length(spc)){
  subset<-anpp_precip%>%
    filter(spc_trt==spc[i])
 test.lm<-lm(anpp~precip_mm, data=subset)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
                      treatment=unique(subset$treatment), 
                      plot_mani=unique(subset$plot_mani), 
                      est=summary(test.lm)$coef["precip_mm", c("Estimate")], 
                      st.er=summary(test.lm)$coef["precip_mm", c("Std. Error")], 
                      p.val=summary(test.lm)$coef["precip_mm","Pr(>|t|)"])
  lm.slopes<-rbind(lm.slopes, output.lm)
}

##test for sig diff between trt-control slopes
##there are so few differences that not going to pay attention to this.

# spc2<-unique(anpp_precip$site_project_comm)
# test.lm<-data.frame()
# for (i in 1:length(spc2)){
#   subset<-anpp_precip%>%
#     filter(site_project_comm==spc2[i])
#   control<-subset%>%
#     filter(plot_mani==0)
#   treat<-subset%>%
#   filter(plot_mani!=0)
# trt_list<-unique(treat$treatment)
# for (i in 1:length(trt_list)){
#   subset2<-treat%>%
#     filter(treatment==trt_list[i])
#   trt<-trt_list[i]
#   ct<-rbind(subset2, control)
#   ct.lm<-lm(anpp~precip_mm*trt, data=ct)
#   output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
#                         treatment=trt, 
#                         est=summary(ct.lm)$coef["precip_mm:trtT", c("Estimate")],
#                         val=summary(ct.lm)$coef["precip_mm:trtT","Pr(>|t|)"])
#   test.lm<-rbind(test.lm, output.lm)
# }
# }

#graphing this
c.slope<-lm.slopes%>%
  filter(plot_mani==0)%>%
  mutate(c_est=est, c_se=st.er)%>%
  select(site_project_comm, c_est, c_se)

t.slope<-lm.slopes%>%
  filter(plot_mani!=0)%>%
  select(-p.val)

slopes_tograph1<-merge(c.slope, t.slope, by="site_project_comm")
slopes_tograph2<-merge(slopes_tograph1, trtint, by=c("site_project_comm","treatment"))
slopes_tograph<-merge(slopes_tograph2, site_info, by="site_project_comm")%>%
  mutate(diff=est-c_est)%>%
  separate(site_project_comm, into=c("site_code","project_name","community_type"), sep="_", remove=F)

###stats
#regression of map with diff
summary(MAP_diff <- lm(diff ~ MAP,  data = slopes_tograph))
#yes sig effect. p = 0.0497
summary(lm(diff ~ MAP,  data = subset(slopes_tograph, trt_type7 == "Nitrogen")))
summary(lm(diff ~ MAP,  data = subset(slopes_tograph, trt_type7 == "Water")))
summary(lm(diff ~ MAP,  data = subset(slopes_tograph, trt_type7 == "Multiple Nutrients")))

#try without MAERC
MAP_diff <- lm(diff ~ MAP,  data = subset(slopes_tograph, site_code!="maerc"))
summary(MAP_diff)
#yes, still sig. p = 0.049


#overall ttest
t.test(slopes_tograph$diff, mu=0)

#do t-test do the slopes differ from zero?
irr<-subset(slopes_tograph, trt_type5=="Water (W)")
t.test(irr$diff, mu=0)

nit<-subset(slopes_tograph, trt_type5=="Nitrogen (N)")
t.test(nit$diff, mu=0)

nuts<-subset(slopes_tograph, trt_type5=="Multiple Nutrients")
t.test(nuts$diff, mu=0)

slopes_bar_trt<-slopes_tograph%>%
  group_by(trt_type6)%>%
  summarise(mdiff=mean(diff),
            ndiff=length(diff),
            sddiff=sd(diff))%>%
  mutate(sediff=sddiff/sqrt(ndiff))%>%
  filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")

slopes_bar_overall<-slopes_tograph%>%
  summarise(mdiff=mean(diff),
            ndiff=length(diff),
            sddiff=sd(diff))%>%
  mutate(sediff=sddiff/sqrt(ndiff))%>%
  mutate(trt_type6="All Treatments")

slopes_bar_overall_abs<-slopes_tograph%>%
  summarise(mdiff=mean(abs(diff)),
            ndiff=length(diff),
            sddiff=sd(abs(diff)))%>%
  mutate(sediff=sddiff/sqrt(ndiff))%>%
  mutate(trt_type6="Abs(All Trts)")

slopes_bar<-rbind(slopes_bar_overall, slopes_bar_trt)
#graphing diff

map<-
ggplot(data=slopes_tograph, aes(x=MAP, y=diff, color = trt_type7))+
  scale_color_manual(name = "GCD Treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green3","darkgray","blue"))+
  geom_point(size=3)+
  geom_smooth(data=subset(slopes_tograph, trt_type6 =="Multiple Nutrients"), method="lm", se=F, color="orange", size = 1)+
  ylab("Difference in Slopes")+
  xlab("Site MAP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", x=275, y=1.6, label="B", size=8)
# 
# control<-
# ggplot(data=slopes_tograph, aes(x=MAP, y=c_est))+
#   geom_point(size=3)+
#   ylab("Slopes for Controls")+
#   xlab("Site MAP")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bar<-
ggplot(data=slopes_bar, aes(x=trt_type6, y=mdiff, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mdiff-sediff, ymax=mdiff+sediff),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Difference in Slopes")+
  #theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_fill_manual(values=c("black", "orange","green3","blue"))+
  xlab("GCD Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.245, label="*", size=8)+
  geom_text(x=2, y=0.24, label="*", size=8)+
  geom_text(x=0.6, y=0.24, label="A", size=8)

grid.arrange(bar, map, ncol=2)

# overall PC anpp ---------------------------------------------------------
all_anpp_dat_mean<-all_anpp_dat%>%
  group_by(site_project_comm, site_code, project_name, community_type, calendar_year, treatment_year, treatment, plot_mani)%>%
  summarize(anpp=mean(anpp))%>%
  ungroup()

controls<-all_anpp_dat_mean%>%
  filter(plot_mani==0)%>%
  mutate(c_anpp=anpp)%>%
  select(site_project_comm, c_anpp, calendar_year, treatment_year)
treatment<-all_anpp_dat_mean%>%
  filter(plot_mani!=0)

fig1<-merge(controls, treatment, by=c("site_project_comm","calendar_year","treatment_year"))%>%
  mutate(anpp_PC=(anpp-c_anpp)/c_anpp)
fig<-merge(fig1, trtint, by=c("site_project_comm",'treatment'))%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))

ggplot(data=fig, aes(x=treatment_year, y=anpp_PC))+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=0.5, color="gray", aes(group=id), se=F)+
  geom_point(size=0.05, color="gray")+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=1, color="black", se=F)+
  geom_hline(yintercept=0)+
  xlab("Treatment Year")+
  ylab("Proportaional Change in ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt_type5, ncol=5, scales="free")



# Biodiversity effects ----------------------------------------------------
#plot level temporal anpp
anpp_temp_cv_plot<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_cv=(anpp_temp_sd/anpp_temp_mean)*100)%>%
  mutate(site_project_comm = paste(site_code, project_name, community_type, sep="_"))

anpp_spc<-all_anpp_dat%>%
  select(site_project_comm)%>%
  unique()

#read in community data
community<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  right_join(anpp_spc)

#get richness for each plot
spc<-unique(community$site_project_comm)
rich_even<-data.frame()

for (i in 1:length(spc)){
  subset<-community%>%
    filter(site_project_comm==spc[i])
  
  out<-community_structure(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id')
  out$site_project_comm<-spc[i]
  
  rich_even<-rbind(rich_even, out)
}

ave_rich<-rich_even%>%
  group_by(site_project_comm, plot_id)%>%
  summarize(richness = mean(richness))%>%
  left_join(anpp_temp_cv_plot)

#1) recreate figures from Yann 2014 Nature paper of NutNet data - not going to present this, there is nothing here.

# controls<-
# ggplot(data=subset(ave_rich, plot_mani==0), aes(x = richness, y = anpp_temp_cv, color = site_project_comm))+
#   geom_point()+
#   xlab("Plot Richness")+
#   ylab("Temporal CV of ANPP")+
#   ggtitle("Control Plots")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
# treated<-
#   ggplot(data=subset(ave_rich, plot_mani!=0), aes(x = richness, y = anpp_temp_cv, color = site_project_comm))+
#   geom_point()+
#   xlab("Plot Richness")+
#   ylab("Temporal CV of ANPP")+
#   ggtitle("Treated Plots")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
# 
# grid.arrange(controls, treated, ncol=2)
# 
#2) recreate Yann 2015 Science paper results
control_rich_stabiltiy<-ave_rich%>%
  filter(plot_mani==0)%>%
  group_by(site_project_comm)%>%
  summarize(cont_rich = mean(richness),
            cont_tempcv = mean(anpp_temp_cv),
            cont_temp_mean = mean(anpp_temp_mean),
            cont_temp_sd = mean(anpp_temp_sd))

stability_logrr<-ave_rich%>%
  filter(plot_mani > 0)%>%
  left_join(control_rich_stabiltiy)%>%
  mutate(log_stability = log(anpp_temp_mean/cont_temp_mean) - log(anpp_temp_sd/cont_temp_sd),
         log_rich = log(richness/cont_rich))%>%
  left_join(trtint)

ggplot(data = stability_logrr, aes(x = log_rich, y = log_stability, color = trt_type7))+
  geom_point()+
  scale_color_manual(name = "GCD Treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  geom_point(size=2)+
  xlab('Change in Richness')+
  ylab('Change in Temporal Stabilty')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)

model_all<-lmer(log_stability ~ log_rich+
                        (0+log_rich | site_code /project_name/community_type/trt_type6),
                      data = stability_logrr)
summary(model_all)

#seperate GCD
model_mult<-lmer(log_stability ~ log_rich+
                  (0+log_rich | site_code /project_name/community_type),
                data = subset(stability_logrr, trt_type7 == "Multiple Nutrients"))
summary(model_mult)
model_water<-lmer(log_stability ~ log_rich+
                   (0+log_rich | site_code /project_name/community_type),
                 data = subset(stability_logrr, trt_type7 == "Water"))
summary(model_water)

model_nit<-lmer(log_stability ~ log_rich+
                   (0+log_rich | site_code /project_name/community_type),
                 data = subset(stability_logrr, trt_type7 == "Nitrogen"))
summary(model_nit)


##doing this for each site with 2 or more experiments
#cdr
model_cdr<-lmer(log_stability ~ log_rich+
                  (0 + log_rich | project_name /community_type/trt_type6),
                data =subset(stability_logrr, site_code == "CDR"))
summary(model_cdr)

#recreating CDR
model_cdr<-lmer(log_stability ~ log_rich*trt_type6+
                  (0 + log_rich|project_name),
                data =subset(stability_logrr, site_code == "CDR"))
summary(model_cdr)
Anova(model_cdr) #yes 

ggplot(data=subset(stability_logrr, site_code == "CDR"), aes(x = log_rich, y = log_stability))+
  geom_point()

#sev
model_sev<-lmer(log_stability ~ log_rich+
                  (log_rich | project_name/community_type/trt_type6),
                data =subset(stability_logrr, site_code == "SEV"))
summary(model_sev)


#knz
model_knz<-lmer(log_stability ~ log_rich+
                  (log_rich | project_name/community_type/trt_type6),
                data =subset(stability_logrr, site_code == "KNZ"))
summary(model_knz)

#recreating konza - yes I generally get the same pattern Mendy is publishing
summary(lm(log_stability ~ log_rich,
        data =subset(stability_logrr, site_code == "KNZ")))
ggplot(data=subset(stability_logrr, site_code == "KNZ"), aes(x = log_rich, y = log_stability))+
  geom_point()

#serc
model_serc<-lmer(log_stability ~ log_rich+
                  (log_rich | project_name / community_type/trt_type6),
                data =subset(stability_logrr, site_code == "SERC"))
summary(model_serc)

ggplot(data = subset(stability_logrr, site_code == "CDR"), aes(x = log_rich, y = log_stability, color = trt_type6))+
  geom_point()+
  geom_smooth(aes(group = trt_type6), method = 'lm', se = F)+
  geom_smooth(method = 'lm', se = F, color = 'black', size = 2)

ggplot(data = subset(stability_logrr, site_code == "KNZ"), aes(x = log_rich, y = log_stability, color = trt_type6))+
  geom_point()+
  geom_smooth(aes(group = trt_type6), method = 'lm', se = F)+
  geom_smooth(method = 'lm', se = F, color = 'black', size = 2)

##3 doing this in a way that makes more sense to me. Instead of doing this for each plot, I would just have one dot for each treatment in each experiment, like all the other figures are done.


###doing for composition change - Not going to present this, I think it does not enhance the paper at all. 
#get composition change
# spc<-unique(community$site_project_comm)
# delta_comp<-data.frame()
# 
# for (i in 1:length(spc)){
#   subset<-community%>%
#     filter(site_project_comm==spc[i])
#   
#   out<-multivariate_change(subset, time.var = 'calendar_year', abundance.var = 'relcov', replicate.var = 'plot_id', treatment = 'treatment', species.var = 'genus_species')
#   out$site_project_comm<-spc[i]
#   
#   delta_comp<-rbind(delta_comp, out)
# }
# ave_comp<-delta_comp%>%
#   group_by(site_project_comm, treatment)%>%
#   summarize(comp_change = mean(composition_change))%>%
#   separate(site_project_comm, into=c("site_code", 'project_name', 'community_type'), sep = "_", remove = F)%>%
#   left_join(anpp_temp_cv)
# 
# control_comp_stabiltiy<-ave_comp%>%
#   filter(plot_mani==0)%>%
#   mutate(cont_comp = comp_change,
#             cont_tempcv = anpp_temp_cv)%>%
#   ungroup()%>%
#   select(-comp_change, -anpp_temp_cv, -treatment, -plot_mani)
# 
# comp_stability_logrr<-ave_comp%>%
#   filter(plot_mani > 0)%>%
#   left_join(control_comp_stabiltiy)%>%
#   mutate(log_stability = log(anpp_temp_cv/cont_tempcv),
#          log_comp = log(comp_change/cont_comp))%>%
#   left_join(trtint)
# 
# 
# ggplot(data = comp_stability_logrr, aes(x = log_comp, y = log_stability, color = trt_type7))+
#   geom_point()+
#   geom_smooth(method = 'lm', se = F, color = 'black', size = 2)+
#   scale_color_manual(name = "GCD treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
#   geom_point(size=2)+
#   geom_smooth(data=subset(comp_stability_logrr, trt_type7 =="Nitrogen"), method="lm", se=F, color="green3", size = 1)+
#   geom_smooth(data=subset(comp_stability_logrr, trt_type7 =="Multiple Nutrients"), method="lm", se=F, color="orange", size = 1)+
#   geom_smooth(data=subset(comp_stability_logrr, trt_type7 =="Water"), method="lm", se=F, color="blue", size = 1)+
#   xlab('Change in Composition')+
#   ylab('Change in Temporal Stabilty')+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##looking into collinearlity in the data
colin_plot<-stability_logrr%>%
  left_join(site_info)

colin_trt<-comp_stability_logrr%>%
  left_join(site_info)
lg

# figure for SEM paper ----------------------------------------------------

####SEM analsyis
sem<-read.csv('SEM_allyr.csv')

sem2<-merge(sem, trtint, by=c("site_project_comm","treatment"))%>%
  mutate(id=paste(site_project_comm, treatment, by="_"))%>%
  filter(trt_type!="CO2"&trt_type!="N+CO2"&trt_type!="N+P+K+irr"&trt_type!="P+K"&trt_type!="drought")

theme_set(theme_bw(12))
ggplot(data=sem2, aes(x=treatment_year, y=anpp_PC, group=trt_type5))+
  #geom_smooth(method="loess", aes(color=trt_type))+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), aes(group=trt_type, color=trt_type), size=1)+
  geom_point(size=0.1, aes(color=trt_type))+
  geom_hline(yintercept=0)+
  facet_wrap(~trt_type, ncol=4, scales="free")


#overall effect of vari --------------------------------------------------
cont_temp<-anpp_temp_cv%>%
  filter(plot_mani==0)%>%
  mutate(cont_temp_cv=anpp_temp_cv)%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, cont_temp_cv)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(-treatment)

trt_temp<-anpp_temp_cv%>%
  filter(plot_mani>0)

tograph1_temp<-merge(cont_temp, trt_temp, by=c("site_code", 'project_name',"community_type"))%>%
  mutate(id=paste(site_code, project_name, community_type, sep="_"))

tograph_temp<-merge(tograph1_temp, trtint, by=c("site_project_comm","treatment"))

##spatial

cont_spat<-lastyr%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_cv=anpp_sp_cv)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_cv, calendar_year)%>%
  select(-treatment)

trt_spat<-lastyr%>%
  filter(plot_mani>0)

tograph1_spat<-merge(cont_spat, trt_spat, by=c("site_project_comm","calendar_year"))

tograph_spat<-merge(tograph1_spat, trtint, by=c("site_project_comm","treatment"))

###other year for spatial
cont_spatother<-other_year%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_cv=anpp_sp_cv)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_cv, calendar_year)%>%
  select(-treatment)

trt_spatother<-other_year%>%
  filter(plot_mani>0)
tograph1_otherspat<-merge(cont_spatother, trt_spatother, by=c("site_project_comm","calendar_year"))

tograph_ottherspat<-merge(tograph1_otherspat, trtint, by=c("site_project_comm","treatment"))

####just overall what are the effects of the treatments on spatial and temporal heterogeneity?
##testing for differences
spat_bar<-tograph_spat%>%
  mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)

t.test(abs(spat_bar$PC), mu=0)
t.test(spat_bar$PC, mu=0)

temp_bar<-tograph_temp%>%
  mutate(PC=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100)

t.test(abs(temp_bar$PC), mu=0)
t.test(temp_bar$PC, mu=0)

##other year
other_bar<-tograph_ottherspat%>%
  mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)


t.test(abs(other_bar$PC), mu=0)

#do t-test do the PC differ from zero?
#spatial
irr<-subset(spat_bar, trt_type6=="Water")
t.test(irr$PC, mu=0)
nit<-subset(spat_bar, trt_type6=="Nitrogen")
t.test(nit$PC, mu=0)
nuts<-subset(spat_bar, trt_type6=="Multiple Nutrients")
t.test(nuts$PC, mu=0)


#temporal model
irr<-subset(temp_bar, trt_type6=="Water")
t.test(irr$PC, mu=0)
nit<-subset(temp_bar, trt_type6=="Nitrogen")
t.test(nit$PC, mu=0)
nuts<-subset(temp_bar, trt_type6=="Multiple Nutrients")
t.test(nuts$PC, mu=0)


tograph_spat_bar_trt<-tograph_spat%>%
  mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)%>%
  group_by(trt_type6)%>%
  summarize(P.C=mean(PC),
            sdd=sd(PC),
            num=length(PC))%>%
  mutate(se=sdd/sqrt(num))%>%
  filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")

tograph_spat_bar_overall<-tograph_spat%>%
  mutate(PC=((anpp_sp_cv-cont_sp_cv)/cont_sp_cv)*100)%>%
  summarize(P.C=mean(PC),
            sdd=sd(PC),
            num=length(PC))%>%
  mutate(se=sdd/sqrt(num))%>%
  mutate(trt_type6="All Treatments")


tograph_spat_bar<-rbind(tograph_spat_bar_overall, tograph_spat_bar_trt)

tograph_temp_trt<-tograph_temp%>%
  mutate(PC=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100)%>%
  group_by(trt_type6)%>%
  summarize(P.C=mean(PC),
            sdd=sd(PC),
            num=length(PC))%>%
  mutate(se=sdd/sqrt(num))%>%
  filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")

tograph_temp_bar_overall<-tograph_temp%>%
  mutate(PC=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100)%>%
  summarize(P.C=mean(PC),
            sdd=sd(PC),
            num=length(PC))%>%
  mutate(se=sdd/sqrt(num))%>%
  mutate(trt_type6="All Treatments")

tograph_temp_bar<-rbind(tograph_temp_trt, tograph_temp_bar_overall)

#graphing this
temp_pc<-
  ggplot(data=tograph_temp_bar, aes(x=trt_type6, y=P.C, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=P.C-se, ymax=P.C+se),position= position_dodge(0.9), width=0.2)+
  ylab("Percent Change of Temporal Variability")+
  scale_fill_manual(values=c("black", "lightgray","gray","darkgray"))+
  xlab("Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  ggtitle("Temporal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")
 
spat_pc<-
  ggplot(data=tograph_spat_bar, aes(x=trt_type6, y=P.C, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=P.C-se, ymax=P.C+se),position= position_dodge(0.9), width=0.2)+
  ylab("Percent Change of Spatial Variability")+
  scale_fill_manual(values=c("black", "lightgray","gray","darkgray"))+
  xlab("Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Spatial")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=40, label="*", size=8)+
  geom_text(x=2, y=60, label="*", size=8)+
  geom_text(x=3, y=85, label="*", size=8)

grid.arrange(temp_pc, spat_pc, ncol=2)



