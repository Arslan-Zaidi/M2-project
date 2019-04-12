#Author AA.Zaidi

#script pulls out putative De Novo somatic heteroplasmies

library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

#read conservative table of heteroplasmies and their counts in all family members
allfam<-fread("files/Analysis/Fixing_heteroplasmy_table/hqcounts_cleaned_conservative_09272018.txt",sep="\t",header=T)

#find potential de novo somatic mutations

#to do this, we first identify sites which are heteroplasmic in only sample of the entire family

#tabulate number of samples that are heteroplasmic in the family
hq.1t.fam<-allfam%>%
  group_by(FID,position)%>%
  summarize(nhets=length(which(maf>0.01)))%>%
  filter(.,nhets==1)%>%
  mutate(fam_het_id=paste(FID,position,sep="_"))

#153 cases

#who are these people
ind.1t<-allfam%>%
  filter(fam_het_id%in%hq.1t.fam$fam_het_id)%>%
  filter(maf>0.01)%>%
  mutate(individual_het_id=paste(individual_id,position,sep="_"))

#make table of frequencies from all family members for these sites
allfam.denovo<-allfam%>%
  filter(fam_het_id%in%hq.1t.fam$fam_het_id)%>%
  mutate(individual_het_id=paste(individual_id,position,sep="_"))

#check if the same site is heteroplasmic in the other tissue
hq.denovo<-allfam.denovo%>%
  filter(fam_het_id%in%allfam.denovo$fam_het_id)%>%
  group_by(FID,fam_het_id,position)%>%
  summarize(nhets.2p=length(which(maf>0.002)))%>%
  filter(nhets.2p==1)

#57 potentially somatic denovo mutations

# #isolate the individual's frequency data
ind.denovo<-allfam.denovo%>%
  filter(fam_het_id%in%hq.denovo$fam_het_id)%>%
  filter(maf>0.01)

#add individuals with no mutations whatsoever - which are also not present in the family
#load complete family file
famfile<-fread("files/Analysis/Fixing_heteroplasmy_table/famfile_cleared_conservative_09272018.txt",header=T,sep="\t")

"%ni%"<-Negate("%in%")

# #assign 0 heteroplasmies to individuals with no denovo mutations
# hq.nodenovo<-famfile%>%
#   filter(individual_id%ni%hq.ndenovo.per.ind$individual_id)%>%
#   group_by(individual_id,tissue,age_collection)%>%
#   summarize(nhets=0)%>%
#   mutate(age_collection=age_collection/365)


#split by tissue
hq.ndenovo.per.ind.bl<-hq.ndenovo.per.ind%>%
  filter(tissue=="bl")

hq.ndenovo.per.ind.ch<-hq.ndenovo.per.ind%>%
  filter(tissue=="ch")

#assign 0 heteroplasmies to individuals with no denovo mutations
hq.nodenovo.bl<-famfile%>%
  filter(individual_id%ni%hq.ndenovo.per.ind.bl$individual_id)%>%
  select(individual_id,age_collection)%>%
  mutate(tissue="bl")%>%
  group_by(individual_id,tissue)%>%
  summarize(nhets=0,age_collection=mean(age_collection))%>%
  mutate(age_collection=age_collection/365)

#assign 0 heteroplasmies to individuals with no denovo mutations
hq.nodenovo.ch<-famfile%>%
  filter(individual_id%ni%hq.ndenovo.per.ind.ch$individual_id)%>%
  select(individual_id,age_collection)%>%
  mutate(tissue="ch")%>%
  group_by(individual_id,tissue)%>%
  summarize(nhets=0,age_collection=mean(age_collection))%>%
  mutate(age_collection=age_collection/365)

#rbind bl
hq.ndenovo.per.ind2.bl<-rbind(hq.ndenovo.per.ind.bl,hq.nodenovo.bl)

#rbind ch
hq.ndenovo.per.ind2.ch<-rbind(hq.ndenovo.per.ind.ch,hq.nodenovo.ch)

#rerun regression


#now separately for each tissue
#blood
hq.age.c2.bl<-glm(data=hq.ndenovo.per.ind2.bl,
                  nhets~age_collection,
                  family="poisson")
hq.ndenovo.per.ind2.bl$pred.nhets=predict(hq.age.c2.bl,type="response")


#cheek
hq.age.c2.ch<-glm(data=hq.ndenovo.per.ind2.ch,
                  nhets~age_collection,
                  family="poisson")

hq.ndenovo.per.ind2.ch$pred.nhets=predict(hq.age.c2.ch,type="response")


hq.ndenovo.per.ind2<-rbind(hq.ndenovo.per.ind2.bl,
                           hq.ndenovo.per.ind2.ch)

bl.beta2<-summary(hq.age.c2.bl)$coefficients[2]
ch.beta2<-summary(hq.age.c2.ch)$coefficients[2]

bl.p2<-summary(hq.age.c2.bl)$coefficients[8]
ch.p2<-summary(hq.age.c2.ch)$coefficients[8]


#make data.frame for glm results
glm.res2<-data.frame(tissue=c("bl","ch"),
                     x=c(50,50),
                     y=c(3.5,3.5),
                     betas=c("Beta= 0.043","Beta= -0.006"),
                     pvalue=c("P= 1.11e-04","P= 0.48"))

#plot relationship between nhets and age at collection
plt.nhets.agec2<-ggplot(hq.ndenovo.per.ind2,aes(age_collection,nhets))+
  geom_point(size=1,alpha=0.4)+
  geom_line(aes(x=age_collection,pred.nhets),color="red")+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        plot.margin = margin(0.2,0.1,0.5,1,"cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank())+
  facet_wrap(~tissue)+
  labs(x="Age at collection of individual",
       title="All ages combined")+
  geom_text(data=glm.res2,aes(x=x,y=y,label=betas))+
  geom_text(data=glm.res2,aes(x=x,y=y-0.5,label=pvalue))


#exclude people < 20 yrs of age
#now separately for each tissue
#blood
hq.ndenovo.per.ind3.bl<-hq.ndenovo.per.ind2.bl%>%
  filter(age_collection>20)
hq.age.c3.bl<-glm(data=hq.ndenovo.per.ind3.bl,
                  nhets~age_collection,
                  family="poisson")
hq.ndenovo.per.ind3.bl$pred.nhets=predict(hq.age.c3.bl,type="response")


#cheek
hq.ndenovo.per.ind3.ch<-hq.ndenovo.per.ind2.ch%>%
  filter(age_collection>20)
hq.age.c3.ch<-glm(data=hq.ndenovo.per.ind3.ch,
                  nhets~age_collection,
                  family="poisson")

hq.ndenovo.per.ind3.ch$pred.nhets=predict(hq.age.c3.ch,type="response")

bl.beta3<-summary(hq.age.c3.bl)$coefficients[2]
ch.beta3<-summary(hq.age.c3.ch)$coefficients[2]

bl.p3<-summary(hq.age.c3.bl)$coefficients[8]
ch.p3<-summary(hq.age.c3.ch)$coefficients[8]


hq.ndenovo.per.ind3<-rbind(hq.ndenovo.per.ind3.bl,
                           hq.ndenovo.per.ind3.ch)

#make data.frame for glm results
glm.res3<-data.frame(tissue=c("bl","ch"),
                    x=c(70,70),
                    y=c(3.5,3.5),
                    betas=c("Beta= 0.04","Beta= -4.2e-04"),
                    pvalue=c("P= 0.029","P= 0.852"))

plt.nhets.agec3<-ggplot(hq.ndenovo.per.ind3,aes(age_collection,nhets))+
  geom_point(size=1,alpha=0.4)+
  geom_line(aes(x=age_collection,pred.nhets),color="red")+
  theme_bw()+
  theme(plot.margin = margin(0,0.1,0.5,1,"cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x=element_blank())+
  facet_wrap(~tissue)+
  labs(x="Age at collection of individual",
       y="No. of de novo mutations in individual's tissues",
       title="Age at collection > 20")+
  geom_text(data=glm.res3,aes(x=x,y=y,label=betas))+
  geom_text(data=glm.res3,aes(x=x,y=y-0.5,label=pvalue))

#check people who are younger than 20
#exclude people > 20 yrs of age
#now separately for each tissue
#blood
hq.ndenovo.per.ind4.bl<-hq.ndenovo.per.ind2.bl%>%
  filter(age_collection<=20)
hq.age.c4.bl<-glm(data=hq.ndenovo.per.ind4.bl,
                  nhets~age_collection,
                  family="poisson")
hq.ndenovo.per.ind4.bl$pred.nhets=predict(hq.age.c4.bl,type="response")


#cheek
hq.ndenovo.per.ind4.ch<-hq.ndenovo.per.ind2.ch%>%
  filter(age_collection<=20)
hq.age.c4.ch<-glm(data=hq.ndenovo.per.ind4.ch,
                  nhets~age_collection,
                  family="poisson")

hq.ndenovo.per.ind4.ch$pred.nhets=predict(hq.age.c4.ch,type="response")

bl.beta4<-summary(hq.age.c4.bl)$coefficients[2]
ch.beta4<-summary(hq.age.c4.ch)$coefficients[2]

bl.p4<-summary(hq.age.c4.bl)$coefficients[8]
ch.p4<-summary(hq.age.c4.ch)$coefficients[8]


hq.ndenovo.per.ind4<-rbind(hq.ndenovo.per.ind4.bl,
                           hq.ndenovo.per.ind4.ch)

#make data.frame for glm results
glm.res4<-data.frame(tissue=c("bl","ch"),
                     x=c(15,15),
                     y=c(3.5,3.5),
                     betas=c("Beta= 0.078","Beta= -0.029"),
                     pvalue=c("P= 0.478","P= 0.478"))

plt.nhets.agec4<-ggplot(hq.ndenovo.per.ind4,aes(age_collection,nhets))+
  geom_point(size=1,alpha=0.4)+
  geom_line(aes(x=age_collection,pred.nhets),color="red")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        plot.margin = margin(0,0.1,0.2,1,"cm"),
        plot.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~tissue)+
  labs(x="Age at collection of individual",
       title="Ages < 20")+
  geom_text(data=glm.res4,aes(x=x,y=y,label=betas))+
  geom_text(data=glm.res4,aes(x=x,y=y-0.5,label=pvalue))
  
plt_all<-plot_grid(plt.nhets.agec2,plt.nhets.agec3,plt.nhets.agec4,
                   ncol=1,
                   align="v",
                   labels=c("A","B","C"),
                   vjust=1.5)

plt_all


ggsave("files/Analysis/Denovo_mutations/plt_somatic_age_04092019.pdf",
       plt_all,
       height=7,
       width=6)
