library(tidyverse)
library(gridExtra)
library(gtable)
library(codyn)
library(lme4)
library(vegan)
library(lmerTest)
library(lmodel2)
library(MVN)

setwd('~/Dropbox/converge_diverge/datasets/LongForm')
setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

theme_set(theme_bw(12)) 

#read in data

# get data ----------------------------------------------------------------


anpp_expInfo<-read.csv("ExperimentInformation_ANPP_Oct2017.csv")%>%
  select(-X)

site_info<-read.csv("SiteExperimentDetails_Dec2016.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, MAP, MAT, rrich, anpp)

anpp<-read.csv("ANPP_Oct2017_2.csv")%>%
  select(-X)%>%
  filter(treatment_year!=0)%>%
  group_by(site_code, project_name, community_type, treatment, plot_id)%>%
  mutate(numyear=length(treatment_year))%>%
  filter(numyear>5)

trtint<-read.csv('treatment interactions_ANPP_datasets_using.csv')%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(site_project_comm, treatment, trt_type7, trt_type5, trt_type6, trt_type)

precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)

precip<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/real_precip_anppSites.csv")%>%
  mutate(calendar_year=year, precip_mm=precip)%>%
  select(-year, -X, -precip)


# clean up anpp data --------------------------------------------------------


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
 

ggplot(data=all_anpp_dat, aes(anpp))+
  geom_histogram()+
  facet_wrap(~site_project_comm, ncol=4, scales="free")


#write.csv(all_anpp_dat, "ANPP_6yrs_Dec2017.csv")
# 
# anpp_trts<-all_anpp_dat%>%
#   select(site_project_comm, treatment)%>%
#   unique


#getting richness_evenness
anpp_spc<-all_anpp_dat%>%
  select(site_project_comm)%>%
  unique()

#read in community data
community<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\SpeciesRelativeAbundance_Oct2017.csv")%>%
  select(-X)%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))%>%
  right_join(anpp_spc)

# calculate temporal cv ---------------------------------------

#Calculate temporal
anpp_temp_stab<-all_anpp_dat%>%
  group_by(site_code, project_name, community_type, treatment,plot_mani, plot_id)%>%
  summarize(anpp_temp_mean=mean(anpp, na.rm=T),
            anpp_temp_sd=sd(anpp, na.rm=T),
            anpp_temp_cv=(anpp_temp_sd/anpp_temp_mean)*100,
            anpp_temp_stab=(anpp_temp_mean/anpp_temp_sd))%>%
  group_by(site_code, project_name, community_type, treatment, plot_mani)%>%
  summarize(anpp_temp_cv=mean(anpp_temp_cv, na.rm=T),
            anpp_temp_mean = mean(anpp_temp_mean, na.rm = T),
            anpp_temp_sd = mean(anpp_temp_sd, na.rm = T),
            anpp_temp_stab= mean(anpp_temp_stab, na.rm=T))

##calculating PC ANPP
meandat<-all_anpp_dat%>%
  group_by(site_project_comm, plot_mani, treatment, treatment_year, calendar_year)%>%
  summarize(manpp=mean(anpp))

mcontrol<-meandat%>%
  filter(plot_mani==0)%>%
  mutate(contanpp=manpp)%>%
  select(-manpp)

mtrt<-meandat%>%
  filter(plot_mani!=0)

PC_anpp<-merge(mtrt, mcontrol, by=c("site_project_comm","treatment_year","calendar_year"))%>%
  mutate(PC=((manpp-contanpp)/contanpp)*100)%>%
  group_by(site_project_comm, treatment.x)%>%
  summarise(mPC=mean(PC))%>%
  mutate(treatment=treatment.x)%>%
  select(-treatment.x)


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
  summarize(varppt=var(precip_mm),
            sdppt=sd(precip_mm))

###calculate community data
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

trt<-community%>%
  select(site_project_comm, plot_id, treatment)%>%
  unique()
plot_mani<-all_anpp_dat%>%
  select(site_project_comm, treatment, plot_mani)%>%
  unique()

ave_rich<-rich_even%>%
  left_join(trt)%>%
  group_by(site_project_comm, plot_id, treatment)%>%
  summarize(richness = mean(richness),
            Evar=mean(Evar, na.rm=T))%>%
  ungroup()%>%
  group_by(site_project_comm, treatment)%>%
  summarize(richness = mean(richness),
            Evar = mean(Evar, na.rm=T))%>%
  right_join(plot_mani)


cont_rich<-ave_rich%>%
  filter(plot_mani==0)%>%
  mutate(cont_rich = richness)%>%
  select(-plot_mani, -treatment, -richness)

PC_rich<-ave_rich%>%
  filter(plot_mani != 0)%>%
  select(-Evar)%>%
  left_join(cont_rich)%>%
  mutate(PC_rich = ((richness-cont_rich)/cont_rich)*100)


# overall effect of vari --------------------------------------------------
cont_temp<-anpp_temp_stab%>%
  filter(plot_mani==0)%>%
  mutate(cont_temp_stab=anpp_temp_stab, 
         cont_temp_mean= anpp_temp_mean,
         cont_temp_sd = anpp_temp_sd,
         cont_temp_cv=anpp_temp_cv)%>%
  ungroup()%>%
  select(site_code, project_name, community_type, treatment, cont_temp_stab, cont_temp_mean,cont_temp_sd, cont_temp_cv)%>%
  mutate(site_project_comm=paste(site_code, project_name,community_type, sep="_"))%>%
  select(-treatment)

trt_temp<-anpp_temp_stab%>%
  filter(plot_mani>0)

CT_comp<-cont_temp%>%
  left_join(trt_temp)%>%
  mutate(id=paste(site_code, project_name, community_type, sep="_"))%>%
  left_join(trtint)%>%
  mutate(PC_stab=(anpp_temp_stab-cont_temp_stab)/cont_temp_stab,
         PC_sd=(anpp_temp_sd-cont_temp_sd)/cont_temp_sd,
         PC_mean=(anpp_temp_mean-cont_temp_mean)/cont_temp_mean,
         PC_CV=(anpp_temp_cv-cont_temp_cv)/cont_temp_cv)


# Q1 how does gcds effect temporal vari? -----------------------

##first overall for PC_CV
with(subset(CT_comp, site_code!="maerc"), t.test(PC_stab, mu=0)) # yes overall sig
t.test(CT_comp$PC_sd, mu=0) #yes overall sig for SD
t.test(CT_comp$PC_mean, mu=0)#yes overall sig for mean
t.test(CT_comp$PC_CV, mu=0)

##num postive or negative PC
sign<-CT_comp%>%
  mutate(pos=ifelse(PC_stab<0, 0,1))

#pos
sum(sign$pos)

48/95 #51% are postive and 49% are negative

##not sig for any
irr<-subset(CT_comp, trt_type6=="Water")
t.test(irr$PC_stab, mu=0)
nit<-subset(CT_comp, trt_type6=="Nitrogen")
t.test(nit$PC_stab, mu=0)
nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
t.test(nuts$PC_stab, mu=0)
## for SD not sig for water or N, sig for mult nuts
irr<-subset(CT_comp, trt_type6=="Water")
t.test(irr$PC_sd, mu=0)
nit<-subset(CT_comp, trt_type6=="Nitrogen")
t.test(nit$PC_sd, mu=0)
nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
t.test(nuts$PC_sd, mu=0)
##for mean sig for all
irr<-subset(CT_comp, trt_type6=="Water")
t.test(irr$PC_mean, mu=0)
nit<-subset(CT_comp, trt_type6=="Nitrogen")
t.test(nit$PC_mean, mu=0)
nuts<-subset(CT_comp, trt_type6=="Multiple Nutrients")
t.test(nuts$PC_mean, mu=0)

PC_bargraph_trt<-CT_comp%>%
  group_by(trt_type6)%>%
  summarize(stab=mean(PC_stab),
            sd_stab=sd(PC_stab),
            sd=mean(PC_sd),
            sd_sd=sd(PC_sd),
            mn=mean(PC_mean),
            sd_mn=sd(PC_mean),
            num=length(PC_CV))%>%
  mutate(se_stab=sd_stab/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  filter(trt_type6=="Nitrogen"|trt_type6=="Multiple Nutrients"|trt_type6=="Water")

CT_comp_trt<-CT_comp%>%
  filter(trt_type7!="Other GCD")
CT_comp_all<-CT_comp%>%
  mutate(trt_type7="All Treatments")
CT_comp_tograph<-rbind(CT_comp_trt, CT_comp_all)

###making boxplots
stab <- ggplot(data = CT_comp_tograph, aes(x = trt_type7, y = PC_stab, color=trt_type7))+
  geom_jitter()+
  geom_boxplot(alpha=.1) +
  xlab("") +
  ylab("Percent Difference\nStability of ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_color_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_hline(yintercept = 0, size = 0.5)+
  geom_text(x=1, y=2.2, label="*", size=8)+
  geom_text(x=0.5, y=1.0, label="C", size=4)

sd <- ggplot(data = CT_comp_tograph, aes(x = trt_type7, y = PC_sd, color=trt_type7))+
  geom_jitter()+
  geom_boxplot(alpha=.1) +
  xlab("") +
  ylab("Percent Difference\nSD of ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_color_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_hline(yintercept = 0, size = 0.5)+
  geom_text(x=1, y=1.7, label="*", size=8)+
  geom_text(x=2, y=1.7, label="*", size=8)+
  geom_text(x=0.5, y=1.7, label="B", size=4)

mean <- ggplot(data = CT_comp_tograph, aes(x = trt_type7, y = PC_mean, color=trt_type7))+
  geom_jitter()+
  geom_boxplot(alpha=.1) +
  xlab("") +
  ylab("Percent Difference\nANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_color_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_hline(yintercept = 0, size = 0.5)+
  geom_text(x=1, y=1.5, label="*", size=8)+
  geom_text(x=2, y=1.5, label="*", size=8)+
  geom_text(x=3, y=.6, label="*", size=8)+
  geom_text(x=4, y=1, label="*", size=8)+
  geom_text(x=0.5, y=1.5, label="A", size=4)

grid.arrange(mean, sd, stab, ncol=1)

##making a bar graph of this
PC_bargraph_all<-CT_comp%>%
  summarize(stab=mean(PC_stab),
            sd_stab=sd(PC_stab),
            sd=mean(PC_sd),
            sd_sd=sd(PC_sd),
            mn=mean(PC_mean),
            sd_mn=sd(PC_mean),
            num=length(PC_CV))%>%
  mutate(se_stab=sd_stab/sqrt(num),
         se_sd=sd_sd/sqrt(num),
         se_mn=sd_mn/sqrt(num))%>%
  mutate(trt_type6="All Treatments")

PC_bargraph<-rbind(PC_bargraph_trt, PC_bargraph_all)


stab_fig<-ggplot(data=PC_bargraph, aes(x=trt_type6, y=stab, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=stab-se_stab, ymax=stab+se_stab),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nCV of ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_fill_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.29, label="*", size=8)+
  geom_text(x=0.6, y=0.3, label="C", size=4)

sd_fig<-ggplot(data=PC_bargraph, aes(x=trt_type6, y=sd, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=sd-se_sd, ymax=sd+se_sd),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nSD of ANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_fill_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+  
  geom_text(x=1, y=0.35, label="*", size=8)+
  geom_text(x=2, y=0.75, label="*", size=8)+
  geom_text(x=0.6, y=0.75, label="B", size=4)+
  scale_y_continuous(limits=c(0, 0.8))

mn_fig<-ggplot(data=PC_bargraph, aes(x=trt_type6, y=mn, fill=trt_type6))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mn-se_mn, ymax=mn+se_mn),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Percent Difference\nANPP")+
  scale_x_discrete(limits = c("All Treatments",'Multiple Nutrients','Nitrogen','Water'),labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("")+
  scale_fill_manual(values=c("orange","green3","blue","black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.35, label="*", size=8)+
  geom_text(x=2, y=0.55, label="*", size=8)+
  geom_text(x=3, y=0.35, label="*", size=8)+
  geom_text(x=4, y=0.48, label="*", size=8)+
  scale_y_continuous(limits=c(0, 0.6))+
  geom_text(x=0.6, y=0.55, label="A", size=4)

grid.arrange(mn_fig, sd_fig, stab_fig, ncol=1)

###what correlates with PC_CV?
#do multiple regression instead of correlations
PC_cor<-CT_comp%>%
  left_join(ave_prod)%>%
  left_join(precip_vari)%>%
  left_join(cont_rich)%>%
  left_join(site_info)%>%
  filter(site_code!="maerc")

tograph_cor<-PC_cor%>%
  select(site_project_comm, treatment, PC_stab, PC_sd, PC_mean, anpp, sdppt, MAP, MAT, cont_rich, Evar)%>%
  gather(parm, value, anpp:Evar)%>%
  gather(vari_metric, vari_value, PC_stab:PC_mean)%>%
  mutate(parm_group=factor(parm, levels = c("cont_rich", "Evar","anpp","MAP","sdppt","MAT")),
         vari_group=factor(vari_metric, levels=c("PC_mean","PC_sd","PC_stab")))

rvalues <- tograph_cor %>% 
  group_by(vari_group, parm_group) %>%
  summarize(r.value = round((cor.test(vari_value, value)$estimate), digits=3),
            p.value = (cor.test(vari_value, value)$p.value))

parameter<-c(
  anpp = "Site ANPP",
  sdppt = "SD of Precip",
  MAP = "MAP",
  MAT = "MAT",
  cont_rich = "Sp Richness",
  Evar = "Evenness"
)

vari<-c(
  PC_stab = "Stability of ANPP",
  PC_mean = "ANPP",
  PC_sd = "SD of ANPP"
)

ggplot(data=tograph_cor, aes(x = value, y = vari_value))+
  geom_point()+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_mean"&parm_group=="Evar"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_mean"&parm_group=="MAT"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_mean"&parm_group=="MAP"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_mean"&parm_group=="anpp"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_sd"&parm_group=="MAP"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_sd"&parm_group=="MAT"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_sd"&parm_group=="sdppt"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_sd"&parm_group=="anpp"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_stab"&parm_group=="anpp"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_stab"&parm_group=="MAP"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_stab"&parm_group=="MAT"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_stab"&parm_group=="sdppt"), method="lm", se=F, color = "black")+
  geom_smooth(data=subset(tograph_cor, vari_group=="PC_stab"&parm_group=="cont_rich"), method="lm", se=F, color = "black")+
  facet_grid(row = vars(vari_group), cols = vars(parm_group), scales="free", labeller=labeller(vari_group = vari, parm_group = parameter))+
  xlab("Value")+
  ylab("Percent Difference")+
  geom_text(data=rvalues, mapping=aes(x=Inf, y = Inf, label = r.value), hjust=1.05, vjust=1.5)


##t-test - do the slopes differ from 1?
#model 2 regression
#first check data is it bivariate normal?
dat<-CT_comp[,c(4,14)]
normal<-mvn(data=log(dat), univariatePlot = "qqplot")#data are bivaiate normal when logged

model2.lm<-lmodel2(log(anpp_temp_stab)~log(cont_temp_stab), range.x = "relative", range.y = "relative", data=dat, nperm=99) #use MA to estimate slope according to package.
#first, I can just use the 97.5% CI interval to say slope does differ from one.
#or I can do a ttest.
slope<-model2.lm$regression.results[2,3]
low<-model2.lm$confidence.intervals[2,4]
high<-model2.lm$confidence.intervals[2,5]
se<-((high-low)/2)/2.24
df<-93
t_value_one <- (slope - 1) / se
2*pt(t_value_one, df=df)
#yes p < 0.029

##doing the same thing for SD not CV
##t-test - do the slopes differ from 1?
#model 2 regression
#first check data is it bivariate normal?
sddat<-CT_comp[,c(6,12)]
normal<-mvn(data=log(sddat), univariatePlot = "qqplot")#data are somewhat bivaiate normal

model2.lm<-lmodel2(log(anpp_temp_sd)~log(cont_temp_sd), range.x = "relative", range.y = "relative", data=sddat, nperm=99) #use MA to estimate slope according to package.
slope<-model2.lm$regression.results[2,3]
low<-model2.lm$confidence.intervals[2,4]
high<-model2.lm$confidence.intervals[2,5]
se<-((high-low)/2)/2.24
df<-93
t_value_one <- (slope - 1) / se
2*pt(t_value_one, df=df)
#yes p < 0.001


# ###looking at three well replicated treatmetns.
# #nitrogen - temporal
subdat<-subset(subset(CT_comp, trt_type6=="Nitrogen"))
dat<-subdat[,c(4,14)]
normal<-mvn(data=dat, univariatePlot = "qqplot")#data are bivaiate normal

model2.lm<-lmodel2(anpp_temp_stab~cont_temp_stab, range.x = "relative", range.y = "relative", data=dat, nperm=99) 
slope<-model2.lm$regression.results[2,3]
low<-model2.lm$confidence.intervals[2,4]
high<-model2.lm$confidence.intervals[2,5]
se<-((high-low)/2)/2.24
df<-9
t_value_one <- (slope - 1) / se
2*pt(t_value_one, df=df)
#not sig.
 
# # nuts temporal
subdat<-subset(subset(CT_comp, trt_type6=="Multiple Nutrients"))
dat<-subdat[,c(4,14)]
normal<-mvn(data=dat, univariatePlot = "qqplot")#data are not and log transfrom doesn't help bivaiate normal

model2.lm<-lmodel2(anpp_temp_stab~cont_temp_stab, range.x = "relative", range.y = "relative", data=dat, nperm=99) 
slope<-model2.lm$regression.results[2,3]
low<-model2.lm$confidence.intervals[2,4]
high<-model2.lm$confidence.intervals[2,5]
se<-((high-low)/2)/2.24
df<-31
t_value_one <- (slope - 1) / se
2*pt(t_value_one, df=df)
#not sig


# #water temporal
subdat<-subset(subset(CT_comp, trt_type6=="Water"))
dat<-subdat[,c(4,14)]
normal<-mvn(data=dat, univariatePlot = "qqplot")#data are not bivaiate normal

model2.lm<-lmodel2(anpp_temp_stab~cont_temp_stab, range.x = "relative", range.y = "relative", data=dat, nperm=99) #use MA to estimate slope according to package.
#first, I can just use the 97.5% CI interval to say slope does differ from one.
#or I can try to do a ttest.
slope<-model2.lm$regression.results[2,3]
low<-model2.lm$confidence.intervals[2,4]
high<-model2.lm$confidence.intervals[2,5]
se<-((high-low)/2)/2.24
df<-5
t_value_one <- (slope - 1) / se
2*pt(t_value_one, df=df, lower=F)
#not sig.

###variance partitioning 
var_temp<- varpart(CT_comp$PC_stab, 
                                ~PC_sd, 
                                ~PC_mean, 
                                data = CT_comp)

### venn diagram plot
plot(var_temp)


###graphing this
theme_set(theme_bw(12))
tograph_color<-CT_comp%>%
  left_join(ave_prod)%>%
  left_join(precip_vari)

##stab figure
dat<-CT_comp[,c(4,14)]
model2.lm<-lmodel2(log(anpp_temp_stab)~log(cont_temp_stab), range.x = "relative", range.y = "relative", data=dat, nperm=99) #use MA 
slopem<-model2.lm$regression.results[2,3]
interceptm<-model2.lm$regression.results[2,2]

#main figure for paper
ggplot(data=CT_comp, aes(x=log(cont_temp_stab), y=log(anpp_temp_stab), color = trt_type7))+
  geom_point(size=3)+
  scale_color_manual(name = "GCD treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  geom_abline(slope=slopem, intercept=interceptm, size=1)+
    ylab("log(Stabilty of Trt Plots)")+
  xlab("log(Stability of Control Plots)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(0,1.9))+
  scale_y_continuous(limits=c(0,1.9))

##sd figure
sddat<-CT_comp[,c(6,12)]
model2.lm<-lmodel2(log(anpp_temp_sd)~log(cont_temp_sd), range.x = "relative", range.y = "relative", data=sddat, nperm=99) #use MA 
slopem<-model2.lm$regression.results[2,3]
interceptm<-model2.lm$regression.results[2,2]

ggplot(data=CT_comp, aes(x=log(cont_temp_sd), y=log(anpp_temp_sd), color=trt_type7))+
  geom_point(size=3)+
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed")+
  scale_color_manual(name = "GCD treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  geom_abline(slope=slopem, intercept=interceptm, size=1)+
  ylab("Log (SD of ANPP Trt Plots)")+
  xlab("Log (SD of ANPP Control Plots)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(3,7))+
  scale_y_continuous(limits=c(3,7))


# Q2 what is the relationship between control CV and PC in ANPP? ---------
C_PC<-PC_anpp%>%
  left_join(cont_temp)%>%
  left_join(trtint)%>%
  left_join(site_info)

#test the relationship between control_temp and PC using model2 regressions
#CV
dat<-C_PC[,c(7,2)]
mvn(data=dat, univariatePlot = "qqplot")
contCV<-lmodel2(mPC~cont_temp_stab, range.x = "relative", range.y = "interval", data=dat, nperm=99)
#Don't use MA bc not the same units. 
cor.test(dat$mPC, dat$cont_temp_stab)#don't use SMA b/c not a significant regression.
hist(dat$mPC)
hist(dat$cont_temp_cv)# no outliers should use RMA
#significat

#SD
sddat<-C_PC[,c(9,2)]
mvn(data=sddat, univariatePlot = "qqplot")
contSD<-lmodel2(mPC~cont_temp_sd, range.x = "relative", range.y = "interval", data=sddat, nperm=99)
#not sig



#looking at three seperate GCDs - NONE are sig.
subdat<-subset(C_PC, trt_type6=="Nitrogen")
dat<-subdat[,c(7,2)]
lmodel2(mPC~cont_temp_stab, range.x = "relative", range.y = "interval", data=dat, nperm=99)

subdat<-subset(C_PC, trt_type6=="Multiple Nutrients")
dat<-subdat[,c(7,2)]
lmodel2(mPC~cont_temp_stab, range.x = "relative", range.y = "interval", data=dat, nperm=99)

subdat<-subset(C_PC, trt_type6=="Water")
dat<-subdat[,c(7,2)]
lmodel2(mPC~cont_temp_stab, range.x = "relative", range.y = "interval", data=dat, nperm=99)

##graphing this

graphQ2<-C_PC%>%
  left_join(ave_prod)%>%
  left_join(precip_vari)

dat<-C_PC[,c(7,2)]
mvn(data=dat, univariatePlot = "qqplot")
contCV<-lmodel2(mPC~cont_temp_stab, range.x = "relative", range.y = "interval", data=dat, nperm=99)
slopem<-contCV$regression.results[4,3]#something is funky with MA USE RMA
interceptm<-contCV$regression.results[4,2]

ggplot(data=graphQ2, aes(x=cont_temp_stab, y=mPC, color = cont_temp_sd, size = cont_temp_mean))+
    geom_point()+
    scale_color_gradient(low = "lightblue", high = "darkred", name = "SD of ANPP")+
    scale_size(name = "Mean ANPP", range = c(1,6))+
  ylab("Percent Difference in ANPP")+
  xlab("Stability of Control Plots")+
  geom_abline(slope=slopem, intercept=interceptm, size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


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

slopes_bar<-rbind(slopes_bar_overall, slopes_bar_trt)
#graphing diff

map<-
ggplot(data=slopes_tograph, aes(x=MAP, y=diff, color = trt_type7))+
  scale_color_manual(name = "GCD Trt", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green3","darkgray","blue"), labels=c("Multiple\nNutrients","Nitrogen","Water","Other GCD"))+
  geom_point(size=3)+
  geom_smooth(method="lm", se=F, color="black", size = 1)+
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
  scale_x_discrete(labels = c("All Trts", "Multiple\n Nutrients", "Nitrogen","Water"))+
  xlab("GCD Treatment")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.245, label="*", size=8)+
  geom_text(x=2, y=0.24, label="*", size=8)+
  geom_text(x=0.6, y=0.24, label="A", size=8)

grid.arrange(bar, map, ncol=2)

bar_poster<-
  ggplot(data=slopes_bar_overall, aes(x=trt_type6, y=mdiff))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=mdiff-sediff, ymax=mdiff+sediff),position= position_dodge(0.9), width=0.2)+
  ylab("")+
  ylab("Difference in Slopes")+
  xlab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 1.5, size = 1)+
  geom_text(x=1, y=0.255, label="*", size=8)


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
  ylab("Percent Difference in ANPP")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~trt_type5, ncol=5, scales="free")



# Biodiversity effects ----------------------------------------------------
#project-treatment level temporal anpp

###relationship between change in richness and change in variabilty
Vari_rich<-PC_rich%>%
  left_join(CT_comp)%>%
  left_join(trtint)

#overall stab
dat<-Vari_rich[,c(7,24)]
mvn(data=dat, univariatePlot = "qqplot")
biod<-lmodel2(PC_stab~PC_rich, range.x = "interval", range.y = "interval", data=dat, nperm=99)
slopem<-biod$regression.results[2,3]
interceptm<-biod$regression.results[2,2]
##overall SD
dat<-Vari_rich[,c(7,23)]
mvn(data=dat, univariatePlot = "qqplot")
biod_sd<-lmodel2(PC_sd~PC_rich, range.x = "interval", range.y = "interval", data=dat, nperm=99)
slopesd<-biod_sd$regression.results[2,3]
interceptsd<-biod_sd$regression.results[2,2]

##treatments seperately
#N not sig
subdat<-subset(subset(Vari_rich, trt_type6=="Nitrogen"))
dat<-subdat[,c(7,24)]
normal<-mvn(data=dat, univariatePlot = "qqplot")#data are bivaiate normal
biod_n<-lmodel2(PC_stab~PC_rich, range.x = "interval", range.y = "interval", data=dat, nperm=99)
#water not sig
subdat<-subset(subset(Vari_rich, trt_type6=="Water"))
dat<-subdat[,c(7,24)]
normal<-mvn(data=dat, univariatePlot = "qqplot")#data are bivaiate normal
biod_w<-lmodel2(PC_stab~PC_rich, range.x = "interval", range.y = "interval", data=dat, nperm=99)
#mult nuts not sig.
subdat<-subset(subset(Vari_rich, trt_type6=="Multiple Nutrients"))
dat<-subdat[,c(7,24)]
normal<-mvn(data=dat, univariatePlot = "qqplot")#data are bivaiate normal
biod_mn<-lmodel2(PC_stab~PC_rich, range.x = "interval", range.y = "interval", data=dat, nperm=99)

#CV pic for paper
ggplot(data = Vari_rich, aes(x = PC_rich, y = PC_stab, color=trt_type7))+
  geom_point()+
  geom_point(size=3)+
  xlab('Percent Difference Richness')+
  ylab('Percent Difference Stability ')+
  scale_color_manual(name = "GCD treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline(slope=slopem, intercept=interceptm, size=1)
#sd for appendix
ggplot(data = Vari_rich, aes(x = PC_rich, y = PC_sd, color=trt_type7))+
  geom_point()+
  geom_point(size=3)+
  xlab('Percent Difference Richness')+
  ylab('Percent Difference SD of ANPP ')+
  scale_color_manual(name = "GCD treatment", breaks = c("Multiple Nutrients","Nitrogen","Water","Other GCD"),values = c("orange", "green2","darkgray","blue"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline(slope=slopesd, intercept=interceptsd, size=1)



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

