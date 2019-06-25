library(tidyverse)
library(gridExtra)
library(codyn)
library(rsq)
library(gtable)
library(grid)
#library(MASS)#loading this package disables select in tidyverse.

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

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
anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
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
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id)%>%
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
cont_temp<-anpp_temp_cv%>%
  filter(plot_mani==0)%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(cont_temp_cv=mean(anpp_temp_cv), 
            cont_temp_mean=mean(anpp_temp_mean),
            cont_temp_sd = mean(anpp_temp_sd))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, cont_temp_cv, cont_temp_mean,cont_temp_sd)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))

trt_temp<-anpp_temp_cv%>%
  filter(plot_mani>0)

CT_comp_plot<-cont_temp%>%
  left_join(trt_temp)%>%
  mutate(id=paste(site_code, project_name, community_type, sep="_"))%>%
  left_join(trtint)%>%
  mutate(PD_CV=((anpp_temp_cv-cont_temp_cv)/cont_temp_cv)*100,
         PD_sd=((anpp_temp_sd-cont_temp_sd)/cont_temp_sd)*100,
         PD_mean=((anpp_temp_mean-cont_temp_mean)/cont_temp_mean)*100)

trt_temp_mean<-anpp_temp_cv%>%
  filter(plot_mani>0)%>%
  group_by(site_code, project_name, community_type, treatment)%>%
  summarize(anpp_temp_cv=mean(anpp_temp_cv), 
            anpp_temp_mean=mean(anpp_temp_mean),
            anpp_temp_sd = mean(anpp_temp_sd))%>%
  ungroup()%>%
  select(site_code, project_name, community_type, anpp_temp_cv, anpp_temp_mean,anpp_temp_sd, treatment)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))


CT_comp<-cont_temp%>%
  left_join(trt_temp_mean)%>%
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
community<-read.csv("SpeciesRelativeAbundance_Oct2017.csv")%>%
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

# Analysis 1. PD diff from zero for each treatment? -----------------------

anpp_temp_cv2<-anpp_temp_cv%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

spc<-unique(anpp_temp_cv2$site_project_comm)
ttest_out<-data.frame()

for (i in 1:length(spc)){
  
  subset<-anpp_temp_cv2%>%
    filter(site_project_comm==spc[i])
  
  trts<-unique(subset(subset, plot_mani!=0)$treatment)
  
  sub_controls<-subset(subset, plot_mani==0)%>%
    rename(cont_mean=anpp_temp_mean,
           cont_cv=anpp_temp_cv)%>%
    ungroup()%>%
    select(cont_mean, cont_cv, site_project_comm)
  
  
  for (j in 1:length(trts)){
    
    sub_treat<-subset(subset, treatment==trts[j])%>%
    ungroup()%>%
      select(anpp_temp_mean, anpp_temp_cv, site_project_comm, treatment)
    
    combined<-sub_treat%>%
      bind_cols(sub_controls)
  
    t.mean = round(t.test(combined$cont_mean, combined$anpp_temp_mean)$statistic, digits=3)
    p.mean = t.test(combined$cont_mean, combined$anpp_temp_mean)$p.value
    cont.mean=t.test(combined$cont_mean, combined$anpp_temp_mean)$estimate[1]
    treat.mean=t.test(combined$cont_mean, combined$anpp_temp_mean)$estimate[2]
    t.cv = round(t.test(combined$cont_cv, combined$anpp_temp_cv)$statistic, digits=3)
    p.cv = t.test(combined$cont_cv, combined$anpp_temp_cv)$p.value
    cont.cv=t.test(combined$cont_cv, combined$anpp_temp_cv)$estimate[1]
    treat.cv=t.test(combined$cont_cv, combined$anpp_temp_cv)$estimate[2]
    
    
    out<-data.frame(site_project_comm=spc[i],
                    treatment=trts[j],
                    t_mean=t.mean, 
                    p_mean=p.mean, 
                    c.mean=cont.mean, 
                    t.mean=treat.mean, 
                    t_cv=t.cv, 
                    p_cv=p.cv, 
                    c.cv=cont.cv, 
                    t.cv=treat.cv)
    
    ttest_out<-rbind(ttest_out, out)
    }
  
}

ttest_summary<-ttest_out%>%
  mutate(meandiff=t.mean-c.mean,
         cvdiff=t.cv-c.cv)%>%
  mutate(resp_mean=ifelse(p_mean>0.05, "not sig", ifelse(p_mean<0.05&meandiff<0, "dec", ifelse(p_mean<0.05&meandiff>0, "inc", 999))),
         resp_cv=ifelse(p_cv>0.05, "not sig", ifelse(p_cv<0.05&cvdiff<0, "dec", ifelse(p_cv<0.05&cvdiff>0, "inc", 999))))%>%
  left_join(trtint)

mean.overall<-ttest_summary%>%
  group_by(resp_mean)%>%
  summarize(n=length(resp_mean))%>%
  mutate(response="A) ANPP")%>%
  rename(effect=resp_mean)%>%
  mutate(trt_type7="All Trts")

mean.trt<-ttest_summary%>%
  group_by(trt_type7, resp_mean)%>%
  summarize(n=length(resp_mean))%>%
  mutate(response="A) ANPP")%>%
  rename(effect=resp_mean)

cv.overall<-ttest_summary%>%
  group_by(resp_cv)%>%
  summarize(n=length(resp_cv))%>%
  mutate(response="B) CV of ANPP")%>%
  rename(effect=resp_cv)%>%
  mutate(trt_type7="All Trts")

cv.trt<-ttest_summary%>%
  group_by(trt_type7, resp_cv)%>%
  summarize(n=length(resp_cv))%>%
  mutate(response="B) CV of ANPP")%>%
  rename(effect=resp_cv)

# Making figure 1 ---------------------------------------------------------
vote.fig<-mean.trt%>%
  bind_rows(cv.trt)%>%
  bind_rows(cv.overall)%>%
  bind_rows(mean.overall)%>%
  mutate(prop=ifelse(trt_type7=="Multiple Nutrients", n/33, ifelse(trt_type7=="Nitrogen", n/11, ifelse(trt_type7=="Water", n/7, ifelse(trt_type7=="Other GCD", n/44, ifelse(trt_type7=="All Trts", n/95, 999))))))

vot<-ggplot(data=vote.fig, aes(y=prop, x=trt_type7, fill=effect))+
  geom_bar(stat="identity")+
  coord_flip()+
  facet_wrap(~response, ncol=1)+
  xlab("Treatment")+
  ylab("Proportion of Treatments Different from Control")+
  scale_fill_manual(name="Treatement Response", label=c("Not Sig.", "Increase", "Decrease"), limits=c("not sig", "inc", "dec"), values = c("Gray", "skyblue", "darkblue"))+
  scale_x_discrete(limits=c("Other GCD", "Water", "Nitrogen", "Multiple Nutrients", "All Trts"))+
  geom_vline(xintercept = 4.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"))

# Analysis 2. is PD different from 0 for each treatment overall -----------

##first overall for PD_CV
t.test(CT_comp$PD_CV, mu=0) # overall No, and not for the difference GCDs
t.test(CT_comp$PD_sd, mu=0) #yes overall sig for SD
t.test(CT_comp$PD_mean, mu=0)#yes overall sig for mean

##num postive or negative PD
sign<-CT_comp%>%
  mutate(pos=ifelse(PD_CV<0, 0,1))

#pos
sum(sign$pos)

45/95 #47% are postive and 53% are negative
50/95
##for not sig for all
irr<-subset(CT_comp, trt_type6=="Water")
t.test(irr$PD_CV, mu=0)
nit<-subset(CT_comp, trt_type6=="Nitrogen")
t.test(nit$PD_CV, mu=0)
nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
t.test(nuts$PD_CV, mu=0)
## for SD sig for all
# irr<-subset(CT_comp, trt_type6=="Water")
# t.test(irr$PD_sd, mu=0)
# nit<-subset(CT_comp, trt_type6=="Nitrogen")
# t.test(nit$PD_sd, mu=0)
# nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
# t.test(nuts$PD_sd, mu=0)
##for mean sig for all
irr<-subset(CT_comp, trt_type6=="Water")
t.test(irr$PD_mean, mu=0)
nit<-subset(CT_comp, trt_type6=="Nitrogen")
t.test(nit$PD_mean, mu=0)
nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
t.test(nuts$PD_mean, mu=0)

###is there a relationship with N level?
nquest<-CT_comp%>%
  left_join(Nlevels)%>%
  filter(trt_type7=="Nitrogen")

summary(lm(PD_CV~n, data=nquest))


# Making figure 2 ---------------------------------------------------------

##making a bar graph of this
PD_bargraph_trt<-CT_comp%>%
  group_by(trt_type6)%>%
  summarize(cv=mean(PD_CV),
            sd_cv=sd(PD_CV),
            sd=mean(PD_sd),
            sd_sd=sd(PD_sd),
            mn=mean(PD_mean),
            sd_mn=sd(PD_mean),
            num=length(PD_CV))%>%
  mutate(se_cv=sd_cv/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")

PD_bargraph_all<-CT_comp%>%
  summarize(cv=mean(PD_CV),
            sd_cv=sd(PD_CV),
            sd=mean(PD_sd),
            sd_sd=sd(PD_sd),
            mn=mean(PD_mean),
            sd_mn=sd(PD_mean),
            num=length(PD_CV))%>%
  mutate(se_cv=sd_cv/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  mutate(trt_type6="All Treatments")

PD_bargraph<-rbind(PD_bargraph_trt, PD_bargraph_all)


cv_fig<-ggplot(data=PD_bargraph, aes(x=trt_type6, y=cv, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=cv-se_cv, ymax=cv+se_cv),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nCV of ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_fill_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=0.6, y=12, label="B", size=4)

# sd_fig<-ggplot(data=PD_bargraph, aes(x=trt_type6, y=sd, fill=trt_type6))+
#   geom_bar(position=position_dodge(), stat="identity")+
#   geom_errorbar(aes(ymin=sd-se_sd, ymax=sd+se_sd),position= position_dodge(0.9), width=0.2)+
#   ylab("")+
#   ylab("Percent Difference\nSD of ANPP")+
#   scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
#   xlab("")+
#   scale_fill_manual(values=c("orange","green3","blue","black"))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
#   geom_vline(xintercept = 1.5, size = 1)+  
#   geom_text(x=1, y=35, label="*", size=8)+
#   geom_text(x=2, y=75, label="*", size=8)+
#   scale_y_continuous(limits=c(0, 80))

mn_fig<-ggplot(data=PD_bargraph, aes(x=trt_type6, y=mn, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mn-se_mn, ymax=mn+se_mn),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_fill_manual(name = "GCD Treatment", values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=35, label="*", size=8)+
  geom_text(x=2, y=55, label="*", size=8)+
  geom_text(x=3, y=35, label="*", size=8)+
  geom_text(x=4, y=48, label="*", size=8)+
  scale_y_continuous(limits=c(0, 60))+
  geom_text(x=0.6, y=55, label="A", size=4)

legend=gtable_filter(ggplot_gtable(ggplot_build(mn_fig)), "guide-box") 
grid.draw(legend)

fig1<-
  grid.arrange(arrangeGrob(mn_fig+theme(legend.position="none"),
                           cv_fig+theme(legend.position="none"),
                           ncol=1), legend, 
               widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)




# Analysis 3, appendix Site-level responses ----------------------------------------------------

####bar graph of difference across ecosystems
PD_ecosystems_test<-CT_comp%>%
  left_join(site_char)%>%
  group_by(site_code, MAP)

cor.test(PD_ecosystems_test$cont_temp_cv, PD_ecosystems_test$manpp)

with(subset(PD_ecosystems_test, site_code!="SERC"), plot(manpp, cont_temp_cv))

ggplot(data=PD_ecosystems_test, aes())

PD_ecosystems<-CT_comp%>%
  left_join(site_char)%>%
  group_by(site_code, MAP)%>%
  summarize(cv=mean(PD_CV),
            sd_cv=sd(PD_CV),
            sd=mean(PD_sd),
            sd_sd=sd(PD_sd),
            mn=mean(PD_mean),
            sd_mn=sd(PD_mean),
            num=length(PD_CV))%>%
  mutate(se_cv=sd_cv/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  ungroup%>%
  mutate(site_code2=ifelse(site_code=="maerc","MAERC", as.character(site_code)))

eco_cv<-ggplot(data=PD_ecosystems, aes(x=reorder(site_code2, MAP), y=cv, fill=MAP))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=cv-se_cv, ymax=cv+se_cv),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nCV of ANPP")+
  xlab("Site Code")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


eco_sd<-ggplot(data=PD_ecosystems, aes(x=reorder(site_code2, MAP), y=sd, fill=MAP))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=sd-se_sd, ymax=sd+se_sd),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nSD of ANPP")+
  xlab("Site Code")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

eco_mn<-ggplot(data=PD_ecosystems, aes(x=reorder(site_code2, MAP), y=mn, fill=MAP))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mn-se_mn, ymax=mn+se_mn),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nANPP")+
  xlab("Site Code")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

legend=gtable_filter(ggplot_gtable(ggplot_build(eco_mn)), "guide-box") 
grid.draw(legend)

fig1<-
  grid.arrange(arrangeGrob(eco_mn+theme(legend.position="none"),
                           eco_cv+theme(legend.position="none"),
                           ncol=1), legend, 
               widths=unit.c(unit(1, "npc") - legend$width, legend$width),nrow=1)

grid.arrange(eco_mn, eco_sd, eco_cv)


# Analysis 4 abiotic and biotic drivers of PD ANPP and CV of ANPP --------------------


###what correlates with PD_CV?

PD_cor<-CT_comp%>%
  left_join(site_char)

#shouldn't use site ANPP as predictive, because it alone is correlated with CV of anpp for the treated and control plots and it is in the 
with(PD_cor, plot(manpp, cont_temp_cv))
with(PD_cor, plot(manpp, anpp_temp_cv))
with(PD_cor, cor.test(manpp, cont_temp_cv))
with(PD_cor, cor.test(manpp, anpp_temp_cv))


#how well does MAP correlated with manpp? Very strongly
with(PD_cor, plot(manpp, MAP))
with(PD_cor, cor.test(manpp, MAP))

#library(MASS) # MASS masks select in tidyverse, so only load this when doing mutliple regressions

##how correlated are the predictor variables?
pairs(PD_cor[,c(21:24,26)])

stepAIC(lm(PD_CV~MAP+MAT+sdppt+rrich+Evar, data=PD_cor))
summary(model.cv<-lm(PD_CV~sdppt+Evar+MAT, data=PD_cor))
rsq.partial(model.cv, adj = T)

# stepAIC(lm(PC_sd~MAT+MAP+anpp+sdppt+cont_rich+Evar, data=PC_cor))
# summary(model.sd<-lm(PC_CV~MAP+anpp+sdppt, data=PC_cor))
# rsq.partial(model.sd, adj =T)

stepAIC(lm(PD_mean~MAT+MAP+sdppt+rrich+Evar, data=PD_cor))
summary(model.mn<-lm(PD_CV~MAT+rrich+Evar, data=PD_cor))
rsq.partial(model.mn)


# Making figure 2 ---------------------------------------------------------

tograph_cor<-PD_cor%>%
  select(site_project_comm, treatment,PD_CV, PD_sd, PD_mean, MAP, sdppt, MAT, rrich, Evar)%>%
  gather(parm, value, MAP:Evar)%>%
  gather(vari_metric, vari_value, PD_CV:PD_mean)%>%
  mutate(parm_group=factor(parm, levels = c("rrich", "Evar","MAP","sdppt","MAT")),
         vari_group=factor(vari_metric, levels=c("PD_mean","PD_sd","PD_CV")))

rvalues <- tograph_cor %>% 
  group_by(vari_group, parm_group) %>%
  summarize(r.value = round((cor.test(vari_value, value)$estimate), digits=3),
            p.value = (cor.test(vari_value, value)$p.value))

parameter<-c(
  MAP = "MAP",
  sdppt = "SD of Precip.",
  MAT = "MAT",
  rrich = "Sp Richness",
  Evar = "Evenness"
)

vari<-c(
  PD_CV = "CV of ANPP",
  PD_mean = "ANPP"
)

tograph_cor2<-tograph_cor%>%
  filter(vari_metric!="PD_sd")
rvalues2<-rvalues %>% 
  filter(vari_group!="PD_sd")

ggplot(data=tograph_cor2, aes(x = value, y = vari_value))+
  geom_point()+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_mean"&parm_group=="Evar"), method="lm", se=F, color = "black")+  
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_mean"&parm_group=="rrich"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_mean"&parm_group=="MAT"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_mean"&parm_group=="MAP"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_mean"&parm_group=="anpp"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_CV"&parm_group=="anpp"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_CV"&parm_group=="MAP"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_CV"&parm_group=="MAT"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_CV"&parm_group=="Evar"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PD_CV"&parm_group=="sdppt"), method="lm", se=F, color = "black")+
  facet_grid(row = vars(vari_group), cols = vars(parm_group), scales="free", labeller=labeller(vari_group = vari, parm_group = parameter))+
  xlab("Value")+
  ylab("Percent Difference")+
  geom_text(data=rvalues2, mapping=aes(x=Inf, y = Inf, label = r.value), hjust=1.05, vjust=1.5)


# Analysis 5 appendix are sites more responsvie to GCDs in low anpp years? ---------

###ARE sites more resopnsive in low ANPP years compared with high ANPP years.

pvalues <- PD_anpp_yr %>% 
  group_by(site_project_comm) %>%
  summarize(p.value = round(summary(lm(PD~contanpp))$coef["contanpp","Pr(>|t|)"], digits=3),
            slope = summary(lm(PD~contanpp))$coef["contanpp", c("Estimate")])%>%
  mutate(pval=ifelse(p.value==0, "<0.001", as.numeric(round(p.value, digits=3))))%>%
  left_join(site_info)

summary(lm(slope~MAP, data=pvalues))
with(pvalues, plot(MAP, slope))

PD_anpp_yr2<-PD_anpp_yr%>%
  mutate(spc_order = factor(site_project_comm, levels = c("SEV_Nfert_0",       "SEV_WENNDEx_0","IMGERS_Yu_0","KLU_KGFert_0","DL_NSFC_0","NWT_snow_0","CDR_BioCON_0" ,"CDR_e001_A","CDR_e001_B","CDR_e001_C", "CDR_e001_D", "CDR_e002_A","CDR_e002_B","CDR_e002_C","KNZ_BGP_0","KNZ_IRG_l","KNZ_IRG_u","KNZ_pplots_0","KNZ_RaMPs_0","KNZ_RHPs_0", "KBS_T7_0","SERC_CXN_0", "SERC_TMECE_MX","SERC_TMECE_SC","SERC_TMECE_SP", "maerc_fireplots_0","ANG_watering_0")))

ggplot(data=PD_anpp_yr2, aes(x = contanpp, y = PD))+
  geom_point()+
  theme(legend.position = "none")+
  facet_wrap(~spc_order, scales = "free")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_BioCON_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e002_B"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_BGP_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_IRG_u"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="maerc_fireplots_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e001_A"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e001_B"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="ANG_watering_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KBS_T7_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e001_C"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e001_D"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_IRG_l"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_RaMPs_0"), method="lm", se=F, color = "black")+
  xlab("Control ANPP")+
  ylab("PD of ANPP")

##do by diff and group by CV
pvalues_df <- PD_anpp_yr %>% 
  group_by(site_project_comm) %>%
  summarize(p.value = round(summary(lm(Diff~contanpp))$coef["contanpp","Pr(>|t|)"], digits=3),
            slope = summary(lm(Diff~contanpp))$coef["contanpp", c("Estimate")])%>%
  mutate(pval=ifelse(p.value==0, "<0.001", as.numeric(round(p.value, digits=3))))%>%
  left_join(site_info)

PD_anpp_yr3<-PD_anpp_yr%>%
  mutate(spc_order = factor(site_project_comm, levels = c("KNZ_RaMPs_0", "KNZ_IRG_l", "KNZ_IRG_u", "KLU_KGFert_0", "SERC_TMECE_MX", "KNZ_pplots_0", "DL_NSFC_0", "ANG_watering_0", "SERC_CXN_0", "SERC_TMECE_SC", "CDR_e002_C", "NWT_snow_0", "CDR_e002_A","SERC_TMECE_SP", "CDR_BioCON_0", "IMGERS_Yu_0", "CDR_e001_A", "KBS_T7_0", "CDR_e001_C", "CDR_e001_D", "CDR_e001_B","CDR_e002_B","SEV_Nfert_0", "SEV_WENNDEx_0", "KNZ_BGP_0", "KNZ_RHPs_0","maerc_fireplots_0")))

ggplot(data=PD_anpp_yr3, aes(x = contanpp, y = Diff))+
  geom_point()+
  theme(legend.position = "none")+
  facet_wrap(~spc_order, scales = "free")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_BioCON_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e002_B"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_BGP_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_IRG_u"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="maerc_fireplots_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="CDR_e001_C"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="KNZ_RaMPs_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="SERC_TMECE_SP"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="IMGERS_Yu_0"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(PD_anpp_yr2, spc_order=="SEV_WENNDEx_0"), method="lm", se=F, color = "black")+
  xlab("Control ANPP")+
  ylab("Diff in ANPP (C-T")

#do this for experiments that are 10 years or longer
PD_anpp_yr_10<-PD_anpp_yr%>%
  filter(treatment_year>10)

summary(aov(lm(PD~treatment_year*site_project_comm, data=PD_anpp_yr_10)))#sig negative slope

ggplot(data=PD_anpp_yr_10, aes(x = treatment_year, y =PD))+
  geom_point()+
  theme(legend.position = "none")+
  geom_smooth(method = "lm")+
  facet_wrap(~site_project_comm, scales = "free")

summary(aov(lm(Diff~contanpp*site_project_comm, data=PD_anpp_yr)))# sig negative slope. Say yes overall negative slope. But differs by sites, X% of studies had a negative slope and there was an interaction between sites. Make fig with p.value in box.

summary(test<-(lm(Diff~contanpp, data=subset(PD_anpp_yr, site_project_comm=="ANG_watering_0"))))

summary(aov(lm(PD~contanpp*site_project_comm, data=PD_anpp_yr)))#sig negative slope

summary(aov(lm(Diff~treatment_year*site_project_comm, data=PD_anpp_yr)))# sig p = 0.048 negative slope

summary(aov(lm(PD~treatment_year*site_project_comm, data=PD_anpp_yr)))#sig negative slope

# Analysis 6, sensitivity of anpp to precip --------------------------------------------------

#precip analysis #1986 in CDR has no precip data, this one year is being dropped.

#drop irrigation treatments becuase it is confusing to add water and then test against precip OR maybe just KNZ because not the same amount of water each year.

#look at average change in sensitivity by treatments.

anpp_precip<-merge(all_anpp_dat, precip, by=c("site_code","calendar_year"))%>%
  mutate(trt=ifelse(plot_mani==0,"C","T"))%>%
  group_by(site_project_comm, trt, calendar_year, precip_mm, treatment, plot_mani)%>%
  summarize(anpp=mean(anpp))%>%
  mutate(spc_trt=paste(site_project_comm, treatment, sep="::"))

# #visually inspecting the data
# ggplot(data=anpp_precip, aes(x=precip_mm, y=anpp, group=treatment, color=trt))+
#   geom_point()+
#   geom_smooth(method="lm", se=F)+
#   facet_wrap(~site_project_comm, ncol=8, scales = "free")


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


#graphing this
c.slope<-lm.slopes%>%
  filter(plot_mani==0)%>%
  rename(c_est=est, c_se=st.er)%>%
  select(site_project_comm, c_est, c_se)

t.slope<-lm.slopes%>%
  filter(plot_mani!=0)%>%
  select(-p.val)

slopes_tograph1<-merge(c.slope, t.slope, by="site_project_comm")
slopes_tograph2<-merge(slopes_tograph1, trtint, by=c("site_project_comm","treatment"))
slopes_tograph<-merge(slopes_tograph2, site_info, by="site_project_comm")%>%
  mutate(diff=est-c_est)%>%
  separate(site_project_comm, into=c("site_code","project_name","community_type"), sep="_", remove=F)%>%
  left_join(precip_vari)

###stats
#regression of map with diff
summary(MAP_diff <- lm(diff ~ MAP,  data = subset(slopes_tograph, MAP<800)))#ns
summary(MAP_diff <- lm(diff ~ MAP,  data = subset(slopes_tograph, MAP>800)))#sig, r=0.11

summary(MAP_diff <- lm(diff ~ MAP,  data = slopes_tograph))
#not sig effect. p = 0.0497
#summary(MAP_diff <- lm(diff ~ sdppt,  data = slopes_tograph))
##yes sig effect. p = 0.0005

summary(lm(diff ~ MAP,  data = subset(slopes_tograph, trt_type7 == "Nitrogen"))) #not sig
summary(lm(diff ~ MAP,  data = subset(slopes_tograph, trt_type7 == "Water")))#not sig
summary(lm(diff ~ MAP,  data = subset(slopes_tograph, trt_type7 == "Multiple Nutrients"))) #sig p = 0.001

#try without MAERC
MAP_diff <- lm(diff ~ MAP,  data = subset(slopes_tograph, site_code!="maerc"&trt_type7=="Multiple Nutrients"))
summary(MAP_diff)
#yes, still sig. p = 0.039


#overall ttest
t.test(slopes_tograph$diff, mu=0)

#do t-test do the slopes differ from zero?
irr<-subset(slopes_tograph, trt_type5=="Water (W)")
t.test(irr$diff, mu=0)

nit<-subset(slopes_tograph, trt_type5=="Nitrogen (N)")
t.test(nit$diff, mu=0)

nuts<-subset(slopes_tograph, trt_type5=="Multiple Nutrients")
t.test(nuts$diff, mu=0)


# Making figure 3 ---------------------------------------------------------
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


slopes_bar<-rbind(slopes_bar_overall, slopes_bar_trt)
#graphing diff

map<-
  ggplot(data=slopes_tograph, aes(x=MAP, y=diff, color = trt_type7))+
  scale_color_manual(name = "GCD Treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green3","darkgray","blue"), labels=c("Multiple\nNutrients","Nitrogen","Water","Other GCD"))+
  geom_point(size=3)+
  # geom_smooth(method="lm", se=F, color="black", size = 1)+
  geom_smooth(data=subset(slopes_tograph, trt_type6 =="Multiple Nutrients"), method="lm", se=F, color="orange", size = 1)+
  ylab("Difference in Slopes")+
  xlab("Site MAP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", x=275, y=1.6, label="B", size=8)

sdppt<-
  ggplot(data=slopes_tograph, aes(x=sdppt, y=diff, color = trt_type7))+
  scale_color_manual(name = "GCD Trt", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green3","darkgray","blue"), labels=c("Multiple\nNutrients","Nitrogen","Water","Other GCD"))+
  geom_point(size=3)+
  geom_smooth(method="lm", se=F, color="black", size = 1)+
  geom_smooth(data=subset(slopes_tograph, trt_type6 =="Multiple Nutrients"), method="lm", se=F, color="orange", size = 1)+
  ylab("Difference in Slopes")+
  xlab("Site SD of Precipitation")+
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
  scale_x_discrete(labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("GCD Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.255, label="*", size=8)+
  geom_text(x=2, y=0.24, label="*", size=8)+
  geom_text(x=0.6, y=0.24, label="A", size=8)

grid.arrange(bar, map, ncol=2)

# analysis 7 appendix is sensitivity related to MAP? ----------------------

##does this differ for wet/dry sites
slopes_bar_site<-slopes_tograph%>%
  group_by(site_code, MAP)%>%
  summarise(mdiff=mean(diff),
            ndiff=length(diff),
            sddiff=sd(diff))%>%
  mutate(sediff=sddiff/sqrt(ndiff))

ggplot(data=slopes_bar_site, aes(x=reorder(toupper(site_code), MAP), y=mdiff, fill=MAP))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mdiff-sediff, ymax=mdiff+sediff),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Difference in Slopes")+
  xlab("Site Code")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# appendix analysis PD anpp over time ---------------------------------------------------------
fig<-PD_anpp_yr%>%
  left_join(trtint)%>%
  mutate(id=paste(site_project_comm, treatment, sep="::"))

ggplot(data=fig, aes(x=treatment_year, y=PD))+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=0.5, color="gray", aes(group=id), se=F)+
  geom_point(size=0.05, color="gray")+
  geom_smooth(method="lm", formula = y ~  poly(x, 2), size=1, color="black", se=F)+
  geom_hline(yintercept=0)+
  xlab("Treatment Year")+
  ylab("Percent Difference in ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt_type5, ncol=5, scales="free")