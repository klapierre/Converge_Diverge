library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(gridExtra)
library(gtools)
library(gtable)
library(grid)
library(lmerTest)

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw(12)) 

#read in data

# get data ----------------------------------------------------------------


anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)

site_info<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, MAP, MAT)

anpp<-read.csv("ANPP_Oct2017_2.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

trtint<-read.csv('treatment interactions_ANPP_datasets_using.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type5, trt_type4, trt_type)

#no longer using prism data
# precip<-read.csv('~/Dropbox/converge_diverge/datasets/LongForm/climate/ANPP_PrecipData.csv')
# 
# precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\ANPP_PrecipData.csv")%>%
#   mutate(site_code=ï..site_code)%>%
#   select(-ï..site_code)

precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)


# clean up anp data --------------------------------------------------------


###select the data to use

#for CDR e001/e002 selecting treatments , 6, 8, 9. For BGP dropping mowing treatments 
#dropping outliers as well
# not using maerc fireplots until I can make more sense of the data

dat2<-merge(anpp_expInfo, anpp, by=c("site_code","project_name","community_type","treatment"))%>%
  select(-nutrients, -light, -carbon, -water, -other_manipulation, -max_trt, -public, -factorial, -block)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  mutate(delete=ifelse(site_code=="CDR"&treatment==2|site_code=="CDR"&treatment==3|site_code=="CDR"&treatment==4|site_code=="CDR"&treatment==5|site_code=="CDR"&treatment==7|site_code=="CDR"&treatment=="2_f_u_n"|site_code=="CDR"&treatment=="3_f_u_n"|site_code=="CDR"&treatment=="4_f_u_n"|site_code=="CDR"&treatment=="5_f_u_n"|site_code=="CDR"&treatment=="7_f_u_n"|project_name=="BGP"&treatment=="u_m_c"|project_name=="BGP"&treatment=="u_m_b"|project_name=="BGP"&treatment=="u_m_n"|project_name=="BGP"&treatment=="u_m_p"|project_name=="BGP"&treatment=="b_m_c"|project_name=="BGP"&treatment=="b_m_b"|project_name=="BGP"&treatment=="b_m_n"|project_name=="BGP"&treatment=="b_m_p"|site_code=="CDR"&anpp>3000|project_name=="BGP"&anpp>2240|project_name=="IRG"&anpp>1500|project_name=="RHPs"&calendar_year==2003, 1, 0))%>%
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

###test another year for spatial data
cont_spat<-other_year%>%
  filter(plot_mani==0)%>%
  mutate(cont_sp_cv=anpp_sp_cv)%>%
  ungroup()%>%
  select(site_project_comm, treatment, cont_sp_cv, calendar_year)%>%
  select(-treatment)

trt_spat<-other_year%>%
  filter(plot_mani>0)
tograph1_otherspat<-merge(cont_spat, trt_spat, by=c("site_project_comm","calendar_year"))

tograph_ottherspat<-merge(tograph1_spat, trtint, by=c("site_project_comm","treatment"))

spat.lm<-lm(anpp_sp_cv~cont_sp_cv, data=tograph_ottherspat)
my.slope <- summary(spat.lm)$coef["cont_sp_cv", c("Estimate", "Std. Error")]
my.df <- summary(spat.lm)$df[2]
t_value_one <- (my.slope["Estimate"] - 1) / my.slope["Std. Error"]
2*pt(t_value_one, df=my.df) # two sided test
# yes p < 0.001

#Does the spatail or temporal CV of control explain variation of treatment?
summary(lm(anpp_temp_cv~cont_temp_cv, data=tograph_temp))
# temp<-lmer(anpp_temp_cv ~ cont_temp_cv +
#                       (cont_temp_cv | site_code / project_name / community_type),
#                     data = tograph_temp)# this allows for slopes and intercetps to vary by experiemnt
# summary(temp)
# anova(temp)
#yes it does

summary(lm(anpp_sp_cv~cont_sp_cv, data=tograph_spat))
# spat<-lmer(anpp_sp_cv ~ cont_sp_cv +
#                (cont_sp_cv | site_code / project_name / community_type),
#              data = tograph_spat)
# summary(spat)
# anova(spat)
#yes it does.

###graphing this
temp<-
ggplot(data=tograph_temp, aes(x=cont_temp_cv, y=anpp_temp_cv))+
  geom_point(aes(color=trt_type5), size=2)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="black")+
  ylab("Temporal CV Treatment Plots")+
  xlab("Temporal CV Control Plots")+
  ggtitle("Temporal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"), breaks=c("CO2 (5)","Irrigation (Irg) (7)","Nitrogen (N) (13)","Phosphorus (7)", "Temperature (Temp) (4)", "Non-Resource (N-R) (7)", "N+CO2 (2)","N+Irg (3)","N+Temp (2)",'Irg+Temp (2)',"Multiple Nutrients (35)","N+Irg+Temp (2)","Nutrients+N-R (6)"))

spat<-
ggplot(data=tograph_spat, aes(x=cont_sp_cv, y=anpp_sp_cv))+
  geom_point(aes(color=trt_type5), size=2)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="black")+
  ylab("Spatial CV Treatment Plots")+
  xlab("Spatial CV Control Plots")+
  ggtitle("Spatial")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"), breaks=c("CO2 (5)","Irrigation (Irg) (7)","Nitrogen (N) (13)","Phosphorus (7)", "Temperature (Temp) (4)", "Non-Resource (N-R) (7)", "N+CO2 (2)","N+Irg (3)","N+Temp (2)",'Irg+Temp (2)',"Multiple Nutrients (35)","N+Irg+Temp (2)","Nutrients+N-R (6)"))

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
  geom_line(aes(group=spc_t, color=trt_type5), size=0.5)+
  ylab("Log RR of Spatial CV")+
  geom_abline(slope=0, intercept=0, size=1)+
  xlab("Treatment Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"), breaks=c("CO2 (5)","Irrigation (Irg) (7)","Nitrogen (N) (13)","Phosphorus (7)", "Temperature (Temp) (4)", "Non-Resource (N-R) (7)", "N+CO2 (2)","N+Irg (3)","N+Temp (2)",'Irg+Temp (2)',"Multiple Nutrients (35)","N+Irg+Temp (2)","Nutrients+N-R (6)"))



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
                      (cont_temp_cv | site_code / project_name / community_type),
                    data = tograph_log_temp)
summary(temp_effect)

#no effect of tempral CV on LogRR, no t-value greater than 2.

#test the relationship between control_spat and effect size
spat_effect <- lmer(mlogrr ~ cont_sp_cv +
                     (cont_sp_cv | site_code / project_name / community_type),
                    data = tograph_log_spat)
summary(spat_effect)


##graphing this
temp_rr<-
  ggplot(data=tograph_log_temp, aes(x=cont_temp_cv, y=mlogrr))+
  geom_point(aes(color=trt_type5), size=2)+
  ylab("Log RR")+
  xlab("Temporal CV Control Plots")+
  ggtitle("Temporal")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"), breaks=c("CO2 (5)","Irrigation (Irg) (7)","Nitrogen (N) (13)","Phosphorus (7)", "Temperature (Temp) (4)", "Non-Resource (N-R) (7)", "N+CO2 (2)","N+Irg (3)","N+Temp (2)",'Irg+Temp (2)',"Multiple Nutrients (35)","N+Irg+Temp (2)","Nutrients+N-R (6)"))


spat_rr<-
ggplot(data=tograph_log_spat, aes(x=cont_sp_cv, y=mlogrr))+
  geom_point(aes(color=trt_type5), size=2)+
  ylab("Log RR")+
  xlab("Spatial CV Control Plots")+
  ggtitle("Spatial")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"), breaks=c("CO2 (5)","Irrigation (Irg) (7)","Nitrogen (N) (13)","Phosphorus (7)", "Temperature (Temp) (4)", "Non-Resource (N-R) (7)", "N+CO2 (2)","N+Irg (3)","N+Temp (2)",'Irg+Temp (2)',"Multiple Nutrients (35)","N+Irg+Temp (2)","Nutrients+N-R (6)"))


legend=gtable_filter(ggplot_gtable(ggplot_build(spat_rr)), "guide-box") 
grid.draw(legend)

grid.arrange(arrangeGrob(temp_rr+theme(legend.position="none"),
                         spat_rr+theme(legend.position="none"),
                         ncol=2), legend, 
             widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)


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

spc2<-unique(anpp_precip$site_project_comm)
test.lm<-data.frame()
for (i in 1:length(spc2)){
  subset<-anpp_precip%>%
    filter(site_project_comm==spc2[i])
  control<-subset%>%
    filter(plot_mani==0)
  treat<-subset%>%
  filter(plot_mani!=0)
trt_list<-unique(treat$treatment)
for (i in 1:length(trt_list)){
  subset2<-treat%>%
    filter(treatment==trt_list[i])
  trt<-trt_list[i]
  ct<-rbind(subset2, control)
  ct.lm<-lm(anpp~precip_mm*trt, data=ct)
  output.lm<-data.frame(site_project_comm=unique(subset$site_project_comm), 
                        treatment=trt, 
                        est=summary(ct.lm)$coef["precip_mm:trtT", c("Estimate")],
                        val=summary(ct.lm)$coef["precip_mm:trtT","Pr(>|t|)"])
  test.lm<-rbind(test.lm, output.lm)
}
}

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

#graphing diff

map<-
ggplot(data=slopes_tograph, aes(x=MAP, y=diff))+
  geom_point(aes(color=trt_type5), size=3)+
  ylab("Difference in Slopes")+
  xlab("Site MAP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(name="Treatment", values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"), breaks=c("CO2 (5)","Irrigation (Irg) (7)","Nitrogen (N) (13)","Phosphorus (7)", "Temperature (Temp) (4)", "Non-Resource (N-R) (7)", "N+CO2 (2)","N+Irg (3)","N+Temp (2)",'Irg+Temp (2)',"Multiple Nutrients (35)","N+Irg+Temp (2)","Nutrients+N-R (6)"))

slopes_bar<-slopes_tograph%>%
  group_by(trt_type5)%>%
  summarise(mdiff=mean(diff),
            ndiff=length(diff),
            sddiff=sd(diff))%>%
  mutate(sediff=sddiff/sqrt(ndiff))%>%
  mutate(trt=ifelse(trt_type5=="CO2 (5)","a",ifelse(trt_type5=="Irrigation (Irg) (7)","b",ifelse(trt_type5=="Nitrogen (N) (13)","c", ifelse(trt_type5=="Phosphorus (7)","d",ifelse(trt_type5=="Temperature (Temp) (4)","e",ifelse(trt_type5=="Non-Resource (N-R) (7)","f",ifelse(trt_type5=="N+CO2 (2)","g",ifelse(trt_type5=="N+Irg (3)","h",ifelse(trt_type5=="N+Temp (2)","i", ifelse(trt_type5=="Irg+Temp (2)","j",ifelse(trt_type5=="Multiple Nutrients (35)","k",ifelse(trt_type5=="N+Irg+Temp (2)","l","m")))))))))))))
  

bar<-
ggplot(data=slopes_bar, aes(x=trt, y=mdiff, fill=trt_type5))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mdiff-sediff, ymax=mdiff+sediff),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_fill_manual(values=c("green","purple","blue","darkred","darkgreen","lightblue","darkorange","yellow3","red","black","gray","pink3","orange"))+
  xlab("Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank())+
  theme(legend.position = "none")


  #annotate(geom="text",x=7, y=.2, label="*", size=8)
  
grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.4, height = 0.4, x = .3, y = 0.3) #plot area for the inset map
print(map,vp=v1) 
print(bar,vp=v2)  


###stats
#do t-test do the slopes differ from zero?
summary(MAP_diff <- lmer(diff ~ trt_type5+
                   (1 | site_code/project_name/community_type),
                 data = slopes_tograph))
Anova(MAP_diff)

co2<-subset(slopes_tograph, trt_type5=="CO2 (5)")
t.test(co2$diff, mu=0)

irr<-subset(slopes_tograph, trt_type5=="Irrigation (Irg) (7)")
t.test(irr$diff, mu=0)

nit<-subset(slopes_tograph, trt_type5=="Nitrogen (N) (13)")
t.test(nit$diff, mu=0)

phos<-subset(slopes_tograph, trt_type5=="Phosphorus (7)")
t.test(phos$diff, mu=0)

temp<-subset(slopes_tograph, trt_type5=="Temperature (Temp) (4)")
t.test(temp$diff, mu=0)

other<-subset(slopes_tograph, trt_type5=="Non-Resource (N-R) (7)")
t.test(other$diff, mu=0)

nco2<-subset(slopes_tograph, trt_type5=="N+CO2 (2)")
t.test(nco2$diff, mu=0)

nirg<-subset(slopes_tograph, trt_type5=="N+Irg (3)")
t.test(nirg$diff, mu=0)

Ntemp<-subset(slopes_tograph, trt_type5=="N+Temp (2)")
t.test(Ntemp$diff, mu=0)

irgtemp<-subset(slopes_tograph, trt_type5=="Irg+Temp (2)")
t.test(irgtemp$diff, mu=0)

nuts<-subset(slopes_tograph, trt_type5=="Multiple Nutrients (35)")
t.test(nuts$diff, mu=0)

nirgtemp<-subset(slopes_tograph, trt_type5=="N+Irg+Temp (2)")
t.test(nirgtemp$diff, mu=0)

nutrn<-subset(slopes_tograph, trt_type5=="Nutrients+N-R (6)")
t.test(nutrn$diff, mu=0)


#regression of map with diff !not working!
MAP_diff <- lmer(diff ~ MAP +
                   trt_type5+
                      (MAP | site_code / project_name / community_type),
                    data = slopes_tograph)
summary(MAP_diff)


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
  mutate(anpp_PC=(c_anpp-anpp)/c_anpp)
fig<-merge(fig1, trtint, by=c("site_project_comm",'treatment'))%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))

ggplot(data=fig, aes(x=treatment_year, y=anpp_PC))+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=0.5, color="black", aes(group=id), se=F)+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=2)+
  geom_point(size=0.1)+
  geom_hline(yintercept=0)+
  xlab("Treatment Year")+
  ylab("Proportaional Change in ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt_type5, ncol=5, scales="free")




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


