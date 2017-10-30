library(tidyverse)
library(ggplot2)

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
#read in data

anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)


anpp<-read.csv("ANPP_Oct2017.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

trtint<-read.csv('treatment interactions_09072017.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type)

dat2<-merge(anpp_expInfo, anpp, by=c("site_code","project_name","community_type","treatment"))%>%
  select(-nutrients, -light, -carbon, -water, -other_manipulation, -max_trt, -public, -factorial, -block)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  mutate(delete=ifelse(site_project_comm=="KNZ_IRG_u"&anpp>1600|site_code=="CDR"&anpp>3000|site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==6|site_code=="CDR"&treatment==8|site_code=="CDR"&treatment=="2_f_u_n"|site_code=="CDR"&treatment=="3_f_u_n"|site_code=="CDR"&treatment=="4_f_u_n"|site_code=="CDR"&treatment=="6_f_u_n"|site_code=="CDR"&treatment=="8_f_u_n",1,0))%>%
  filter(delete!=1)

nosev<-dat2%>%
  filter(site_project_comm!="SEV_Nfert_0")

sev<-dat2%>%
  filter(site_project_comm=="SEV_Nfert_0")%>%
  select(-treatment_year)%>%
  mutate(treatment_year= ifelse(calendar_year==2004, 10, ifelse(calendar_year==2005, 11, ifelse(calendar_year==2006, 12, ifelse(calendar_year==2007, 13, ifelse(calendar_year==2008, 14, ifelse(calendar_year==2009, 15, ifelse(calendar_year==2010, 16, ifelse(calendar_year==2011, 17, 18)))))))))

all_anpp_dat<-rbind(sev, nosev)

ggplot(data=all_anpp_dat, aes(anpp))+
  geom_histogram()+
  facet_wrap(~site_project_comm, ncol=4, scales="free")


#write.csv(all_anpp_dat, "ANPP_6yrs_alldata_fixedsev_cdrsubset.csv")

#calculate spatail
anpp_spatial<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment_year, calendar_year, treatment, plot_mani)%>%
  summarize(anpp_sp_sd=sd(anpp, na.rm=T),
            anpp_sp_mean=mean(anpp, na.rm=T),
            anpp_sp_cv=(anpp_sp_sd/anpp_sp_mean)*100)%>%
  select(-anpp_sp_sd, -anpp_sp_mean)%>%
  filter(treatment_year!=0)

#Calculate temporal
anpp_temp_cv<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_cv=(anpp_temp_sd/anpp_temp_mean)*100)%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarize(anpp_temp_cv=mean(anpp_temp_cv, na.rm=T))

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

logRRsp<-logRR<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(logrr=abs(log(manpp/contanpp)))%>%
  group_by(site_project_comm, treatment.x, calendar_year)%>%
  summarise(mlogrr=mean(logrr),
            sdlogrr=sd(logrr))%>%
  mutate(treatment=treatment.x)%>%
  ungroup()%>%
  select(-treatment.x)

###temporal analysis
cont<-anpp_temp_cv%>%
  filter(plot_mani==0)%>%
  mutate(cont_temp_cv=anpp_temp_cv)%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, cont_temp_cv)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(-treatment)

trt<-anpp_temp_cv%>%
  filter(plot_mani>0)

tograph1<-merge(cont, trt, by=c("site_code", 'project_name',"community_type"))%>%
  mutate(id=paste(site_code, project_name, community_type, sep="_"))

tograph<-merge(tograph1, trtint, by=c("site_project_comm","treatment"))

ggplot(data=tograph, aes(x=cont_temp_cv, y=anpp_temp_cv))+
  geom_point(aes(color=trt_type))+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method="lm")

tograph_log1<-merge(logRR, cont, by="site_project_comm")
tograph_log<-merge(tograph_log1, trtint, by=c("site_project_comm","treatment"))

ggplot(data=tograph_log, aes(x=cont_temp_cv, y=mrr))+
  geom_point(aes(color=trt_type))+
  geom_smooth(method='lm')



###spatial analysis
lastyr<-anpp_spatial%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  group_by(site_project_comm)%>%
  mutate(maxyr=max(calendar_year))%>%
  filter(maxyr==calendar_year)

cont<-lastyr%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_cv=anpp_sp_cv)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_cv, calendar_year)%>%
  select(-treatment)

trt<-lastyr%>%
  filter(plot_mani>0)

tograph1<-merge(cont, trt, by=c("site_project_comm","calendar_year"))

tograph<-merge(tograph1, trtint, by=c("site_project_comm","treatment"))

ggplot(data=tograph, aes(x=cont_sp_cv, y=anpp_sp_cv))+
  geom_point(aes(color=trt_type))+
  geom_abline(slope=1, intercept=0)+
  geom_smooth(method="lm")

tograph_log1<-merge(logRRsp, cont, by=c("site_project_comm","calendar_year"))
tograph_log<-merge(tograph_log1, trtint, by=c("site_project_comm","treatment"))

ggplot(data=tograph_log, aes(x=cont_sp_cv, y=mlogrr))+
  geom_point(aes(color=trt_type))+
  geom_smooth(method="lm")



####SEM analsyis
sem<-read.csv('SEM_allyr.csv')

sem2<-merge(sem, trtint, by=c("site_project_comm","treatment"))%>%
  mutate(id=paste(site_project_comm, treatment, by="_"))%>%
  filter(trt_type!="CO2"&trt_type!="N+CO2"&trt_type!="N+P+K+irr"&trt_type!="P+K"&trt_type!="drought")

theme_set(theme_bw(12))
ggplot(data=sem2, aes(x=treatment_year, y=anpp_PC, group=trt_type))+
  #geom_smooth(method="loess", aes(color=trt_type))+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), aes(group=trt_type, color=trt_type), size=1)+
  geom_point(size=0.1, aes(color=trt_type))+
  geom_hline(yintercept=0)+
  facet_wrap(~trt_type, ncol=4, scales="free")


