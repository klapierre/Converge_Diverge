library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(gridExtra)
library(gtools)
library(gtable)
library(grid)

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw(12))

#read in data

# get data ----------------------------------------------------------------


anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)


anpp<-read.csv("ANPP_Oct2017_2.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

trtint<-read.csv('treatment interactions_ANPP_datasets_using.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type, trt_type2, trt_type3)

precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\ANPP_fromPrism.csv")%>%
  mutate(site_code=ï..site_code)%>%
  select(-ï..site_code)

# clean up anp data --------------------------------------------------------


###select the data to use

#for CDR e001/e002 selecting treatments , 6, 8, 9. For BGP dropping mowing treatments 
#dropping outliers as well

dat2<-merge(anpp_expInfo, anpp, by=c("site_code","project_name","community_type","treatment"))%>%
  select(-nutrients, -light, -carbon, -water, -other_manipulation, -max_trt, -public, -factorial, -block)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  mutate(delete=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7|site_code=="CDR"&treatment=="2_f_u_n"|site_code=="CDR"&treatment=="3_f_u_n"|site_code=="CDR"&treatment=="4_f_u_n"|site_code=="CDR"&treatment=="5_f_u_n"|site_code=="CDR"&treatment=="7_f_u_n"|project_name=="BGP"&treatment=="u_m_c"|project_name=="BGP"&treatment=="u_m_b"|project_name=="BGP"&treatment=="u_m_n"|project_name=="BGP"&treatment=="u_m_p"|project_name=="BGP"&treatment=="b_m_c"|project_name=="BGP"&treatment=="b_m_b"|project_name=="BGP"&treatment=="b_m_n"|project_name=="BGP"&treatment=="b_m_p"|site_project_comm=="maerc_fireplots_0"&anpp>3500|site_code=="CDR"&anpp>3000|project_name=="BGP"&anpp>2500|project_name=="IRG"&anpp>1500,1,0))%>%
  filter(delete!=1)

nosev<-dat2%>%
  filter(site_project_comm!="SEV_Nfert_0")

sev<-dat2%>%
  filter(site_project_comm=="SEV_Nfert_0")%>%
  select(-treatment_year)%>%
  mutate(treatment_year= ifelse(calendar_year==2004, 10, ifelse(calendar_year==2005, 11, ifelse(calendar_year==2006, 12, ifelse(calendar_year==2007, 13, ifelse(calendar_year==2008, 14, ifelse(calendar_year==2009, 15, ifelse(calendar_year==2010, 16, ifelse(calendar_year==2011, 17, 18)))))))))

all_anpp_dat<-rbind(sev, nosev)

sites<-all_anpp_dat%>%
  select(site_project_comm, treatment, plot_mani)%>%
  unique()%>%
  filter(plot_mani!=0)

ggplot(data=all_anpp_dat, aes(anpp))+
  geom_histogram()+
  facet_wrap(~site_project_comm, ncol=4, scales="free")


write.csv(all_anpp_dat, "ANPP_6yrs_Oct2017.csv")
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
            anpp_sp_cv=(anpp_sp_sd/anpp_sp_mean)*100)%>%
  select(-anpp_sp_sd, -anpp_sp_mean)%>%
  filter(treatment_year!=0)

lastyr<-anpp_spatial%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  group_by(site_project_comm)%>%
  mutate(maxyr=max(calendar_year))%>%
  filter(maxyr==calendar_year)

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

logRRsp<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(logrr=abs(log(manpp/contanpp)))%>%
  group_by(site_project_comm, treatment.x, calendar_year)%>%
  summarise(mlogrr=mean(logrr),
            sdlogrr=sd(logrr))%>%
  mutate(treatment=treatment.x)%>%
  ungroup()%>%
  select(-treatment.x)



# Q1 how does gcds effect temporal or spatial vari? -----------------------

###temporal analysis
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

##t-test - do the slopes differ from 1?
#Temporal control versus treatment
# do the lines differ from a slope of 1?
temp.lm<-lm(anpp_temp_cv~cont_temp_cv, data=tograph_temp)
my.slope <- summary(temp.lm)$coef["cont_temp_cv", c("Estimate", "Std. Error")]
my.df <- summary(temp.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001

#spatail control versus treatment
# do the lines differ from a slope of 1?
spat.lm<-lm(anpp_sp_cv~cont_sp_cv, data=tograph_spat)
my.slope <- summary(spat.lm)$coef["cont_sp_cv", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001

###graphing this
temp<-
ggplot(data=tograph_temp, aes(x=cont_temp_cv, y=anpp_temp_cv))+
  geom_point(aes(color=trt_type3), size=2)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="black")+
  ylab("Temporal CV Treatment Plots")+
  xlab("Temporal CV Control Plots")+
  ggtitle("Temporal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("lightblue","green","blue","pink3","black","red","gray","purple","darkgray","red"), breaks=c("CO2","Irrigation","Nitrogen","Phosphorus","Other", "2 Resources","Multiple Resources","Resource+Other","Multiple Resources+Other"))

spat<-
ggplot(data=tograph_spat, aes(x=cont_sp_cv, y=anpp_sp_cv))+
  geom_point(aes(color=trt_type3), size=2)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="black")+
  ylab("Spatial CV Treatment Plots")+
  xlab("Spatial CV Control Plots")+
  ggtitle("Spatial")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("lightblue","green","blue","pink3","black","red","gray","purple","darkgray","red"), breaks=c("CO2","Irrigation","Nitrogen","Phosphorus","Other", "2 Resources","Multiple Resources","Resource+Other","Multiple Resources+Other"))

legend=gtable_filter(ggplot_gtable(ggplot_build(spat)), "guide-box") 
grid.draw(legend)

grid.arrange(arrangeGrob(temp+theme(legend.position="none"),
                         spat+theme(legend.position="none"),
                         ncol=2), legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)


###spatial through time.
cont_spat_all<-anpp_spatial%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_cv=anpp_sp_cv)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_cv, calendar_year, treatment_year)%>%
  select(-treatment)

trt_spat_all<-anpp_spatial%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  filter(plot_mani>0)

tograph1_spat_all<-merge(cont_spat_all, trt_spat_all, by=c("site_project_comm","treatment_year","calendar_year"))

tograph_spat_all<-merge(tograph1_spat_all, trtint, by=c("site_project_comm","treatment"))%>%
  mutate(logrr=log(anpp_sp_cv/cont_sp_cv))%>%
  mutate(spc_t=paste(site_project_comm, treatment, sep=""))

ggplot(data=tograph_spat_all, aes(x=treatment_year, y=logrr))+
  geom_line(aes(group=spc_t, color=trt_type3), size=0.5)+
  ylab("Log RR of Spatial CV")+
  geom_abline(slope=0, intercept=0, size=1)+
  xlab("Treatment Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("lightblue","green","blue","pink3","black","red","gray","purple","darkgray","red"), breaks=c("CO2","Irrigation","Nitrogen","Phosphorus","Other", "2 Resources","Multiple Resources","Resource+Other","Multiple Resources+Other"))


# Q2 what is the relationship between control CV and effect size? ---------

##temporal
tograph_log1_temp<-merge(logRR, cont_temp, by="site_project_comm")
tograph_log_temp<-merge(tograph_log1_temp, trtint, by=c("site_project_comm","treatment"))

###spatial analysis
tograph_log1_spat<-merge(logRRsp, cont_spat, by=c("site_project_comm","calendar_year"))
tograph_log_spat<-merge(tograph_log1_spat, trtint, by=c("site_project_comm","treatment"))%>%
  separate(site_project_comm, into=c("site_code","project_name","community_type"), sep="_", remove=F)


#mixed-model
#test the relationship between control_temp and effect size
temp_effect <- lmer(mlogrr ~ cont_temp_cv +
                      trt_type3+
                      (cont_temp_cv | site_code / project_name / community_type),
                    data = tograph_log_temp)
summary(temp_effect)

#no effect of tempral CV on LogRR, no t-value greater than 2.

#test the relationship between control_spat and effect size
spat_effect <- lmer(mlogrr ~ cont_sp_cv +
                      trt_type3+
                      (cont_sp_cv | site_code / project_name / community_type),
                    data = tograph_log_spat)
summary(spat_effect)


##graphing this
temp_rr<-
  ggplot(data=tograph_log_temp, aes(x=cont_temp_cv, y=mlogrr))+
  geom_point(aes(color=trt_type3), size=2)+
  ylab("Log RR")+
  xlab("Temporal CV Control Plots")+
  ggtitle("Temporal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("lightblue","green","blue","pink3","black","red","gray","purple","darkgray","red"), breaks=c("CO2","Irrigation","Nitrogen","Phosphorus","Other", "2 Resources","Multiple Resources","Resource+Other","Multiple Resources+Other"))

spat_rr<-
ggplot(data=tograph_log_spat, aes(x=cont_sp_cv, y=mlogrr))+
  geom_point(aes(color=trt_type3), size=2)+
  ylab("Log RR")+
  xlab("Spatial CV Control Plots")+
  ggtitle("Spatial")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("lightblue","green","blue","pink3","black","red","gray","purple","darkgray","red"), breaks=c("CO2","Irrigation","Nitrogen","Phosphorus","Other", "2 Resources","Multiple Resources","Resource+Other","Multiple Resources+Other"))


legend=gtable_filter(ggplot_gtable(ggplot_build(spat_rr)), "guide-box") 
grid.draw(legend)

grid.arrange(arrangeGrob(temp_rr+theme(legend.position="none"),
                         spat_rr+theme(legend.position="none"),
                         ncol=2), legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)


# precipitation analysis --------------------------------------------------

#precip analysis
#this will drop experiments at sites KLU, DL, IMGERS, a total of 3 experiments because only have data from US sites

anpp_precip<-merge(all_anpp_dat, precip, by=c("site_code","calendar_year"))%>%
  mutate(trt=ifelse(plot_mani==0,"C","T"))%>%
  group_by(site_project_comm, trt, calendar_year, ppt_mm, treatment, plot_mani)%>%
  summarize(anpp=mean(anpp))%>%
  mutate(spc_trt=paste(site_project_comm, treatment, sep="::"))


ggplot(data=anpp_precip, aes(x=ppt_mm, y=anpp, group=treatment, color=trt))+
  geom_point()+
  geom_smooth(method="lm", se=F)+
  facet_wrap(~site_project_comm, ncol=8, scales = "free")

spc<-unique(anpp_precip$spc_trt)

lm.slopes<-data.frame()

for (i in 1:length(spc)){
  subset<-anpp_precip%>%
    filter(spc_trt==spc[i])
  
 test.lm<-lm(anpp~ppt_mm, data=subset)
  
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
                      trt=unique(subset$trt), 
                      treatment=unique(subset$treatment), 
                      plot_mani=unique(subset$plot_mani), 
                      est=summary(test.lm)$coef["ppt_mm", c("Estimate")], 
                      st.er=summary(test.lm)$coef["ppt_mm", c("Std. Error")], 
                      p.val=summary(test.lm)$coef["ppt_mm","Pr(>|t|)"])
  
  lm.slopes<-rbind(lm.slopes, output.lm)
}
  
treat.lm<-data.frame()

for (i in 1:length(spc)){
  subset<-anpp_precip%>%
    filter(site_project_comm==spc[i])%>%
    filter(plot_mani!=0)
  
  control<-subset%>%
    filter(plot_mani==0)
  
  c.lm<-lm(anpp~ppt_mm, data=control)
  
  cont.lm<-data.frame(site_project_comm=unique(control$site_project_comm), 
                      trt=unique(control$trt), 
                      treatment=unique(control$treatment), 
                      plot_mani=unique(control$plot_mani), 
                      est=summary(c.lm)$coef["ppt_mm", c("Estimate")], 
                      st.er=summary(c.lm)$coef["ppt_mm", c("Std. Error")], 
                      p.val=summary(c.lm)$coef["ppt_mm","Pr(>|t|)"])
  
  control.lm<-rbind(control.lm, cont.lm)
}
  treat<-subset%>%
    filter(plot_mani!=0)
  
  trt_list<-unique(treat$treatment)
  
  for (i in 1:length(trt_list)){
    
    subset2<-treat%>%
      filter(treatment==trt_list[i])
    
  ct<-rbind(subset2, control)
  
  ct.lm<-lm(anpp~ppt_mm*trt, data=ct)
     
  }
  
}


# figure for SEM paper ----------------------------------------------------

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


