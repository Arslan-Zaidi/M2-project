#Author AA.Zaidi

#script pulls out putative De Novo germline heteroplasmies

library(data.table)
library(dplyr)
library(ggplot2)

#read conservative table of heteroplasmies and their counts in all family members
allfam<-fread("files/Analysis/Fixing_heteroplasmy_table/hqcounts_cleaned_conservative_09272018.txt",sep="\t",header=T)

#find potential de novo mutations

#to do this, we first identify samples who are heteroplasmic for both tissues but no one in their family shares that heteroplasmy

#tabulate number of samples that are heteroplasmic in the family
hq.2t.fam<-allfam%>%
  group_by(FID,position)%>%
  summarize(nhets=length(which(maf>0.01)))%>%
  filter(.,nhets==2)%>%
  mutate(fam_het_id=paste(FID,position,sep="_"))

#subset cases in which both of these heteroplasmies are present in the same individual
hq.2t.ind<-allfam%>%
  filter(fam_het_id%in%hq.2t.fam$fam_het_id)%>%
  group_by(FID,fam_het_id,fam_cat,individual_id,level,position)%>%
  summarize(nhets=length(which(maf>0.01)))%>%
  filter(nhets==2)

#select cases only where these are present in kids
hq.denovo<-hq.2t.ind%>%
  filter(
    fam_cat%in%c("g1m1c2","tg1m1c1c2g2","g1m1c2m2c0","g1m1c3") & level=="m1"|
    fam_cat%in%c("m1c2","m1c3","m1c4","m1c5") & level%in%c("c1","c2","c3","c4","c5"))%>%
  mutate(individual_het_id=paste(individual_id,position,sep="_"))

#103 cases

#extract data for sites which are denovo from all family members
allfam.denovo<-allfam%>%
  filter(fam_het_id%in%hq.denovo$fam_het_id)%>%
  mutate(individual_het_id=paste(individual_id,position,sep="_"))


#count number of individuals in each family who are not the "proband" and have heteroplasmy at that site less than 0.2%
#only keep instances where the remaining family members have a MAF<0.2%
hq.denovo.red<-allfam.denovo%>%
  group_by(FID,position)%>%
  summarize(nhets.2p=length(setdiff(which(maf>0.002),which(individual_het_id%in%hq.denovo$individual_het_id))))%>%
  mutate(fam_het_id=paste(FID,position,sep="_"))%>%
  filter(nhets.2p==0)

#78 cases

#isolate these cases - frequency info for all family members
allfam.denovo.red<-allfam%>%
  filter(fam_het_id%in%hq.denovo.red$fam_het_id)%>%
  mutate(individual_het_id=paste(individual_id,position,sep="_"))

#output these heteroplasmies for supplement
hq.ind.denovo<-allfam.denovo.red%>%
  filter(maf>0.01)%>%
  select(FID,mother_id,individual_id,level,tissue,position,maf)

fwrite(hq.ind.denovo,"files/Analysis/Denovo_mutations/hq_denovo_mutations_03202019.txt",sep="\t",col.names=T,row.names=F,quote=F)

#count the number of de novo heteroplasmies per individual
hq.ndenovo.per.ind<-allfam.denovo.red%>%
  filter(maf>0.01)%>%
  group_by(FID,level,individual_id,age_collection,age_birth)%>%
  summarise(ndenovo.hets=length(unique(position)))%>%
  group_by(FID,level,individual_id,age_collection,age_birth)%>%
  summarise(nhets=sum(ndenovo.hets))%>%
  mutate(age.c=age_collection/365,age.b=age_birth/365)

"%ni%"<-Negate("%in%")

#now add 0s - the kids who don't have any denovo mutations

#for this readin the whole famfile 
famfile<-fread("files/Analysis/Fixing_heteroplasmy_table/famfile_cleared_conservative_09272018.txt",header=T,sep="\t")

allfam.nodenovo<-famfile[which((famfile$individual_id%ni%hq.ndenovo.per.ind$individual_id &
                   famfile$fam_cat%in%c("g1m2","g1m1c2","g1m1c2m2c0","g1m1c3") & 
                   famfile$level%in%c("m1","m2","c1","c2","c3"))|
                     (famfile$individual_id%ni%hq.ndenovo.per.ind$individual_id &
                        famfile$fam_cat%in%c("tg1m1c3","tg1m1c1c2g2") & 
                        famfile$level%in%c("g1","g2","m1","c1","c2","c3"))|
                     (famfile$individual_id%ni%hq.ndenovo.per.ind$individual_id &
                        famfile$fam_cat%in%c("m1c2","m1c3","m1c4","m1c5") & 
                        famfile$level%in%c("c1","c2","c3","c4","c5"))),]


nodenovo.per.ind<-allfam.nodenovo%>%
  group_by(FID,level,individual_id,age_collection,age_birth)%>%
  summarize(nhets=0)%>%
  mutate(age.c=age_collection/365,age.b=age_birth/365)


hq.ndenovo.per.ind<-rbind(hq.ndenovo.per.ind,nodenovo.per.ind)

hq.ndenovo.per.ind<-hq.ndenovo.per.ind[order(hq.ndenovo.per.ind$individual_id),]

#poisson regression of nhets ~ age at collection - defunct no longer valid hypothesis
# hq.age.c<-glm(data=hq.ndenovo.per.ind,nhets~age.c,family="poisson")
# pred.age.c<-data.frame(age.c=hq.ndenovo.per.ind$age.c,pred.nhets=predict(hq.age.c,type="response"))

#poisson regression of nhets ~ age at birth
pred.age.b<-hq.ndenovo.per.ind%>%
  filter(is.na(age.b)=="FALSE")
hq.age.b<-glm(data=pred.age.b,nhets~age.b,family="poisson")
pred.age.b$pred.nhets=predict(hq.age.b,type="response")

#plot relationship between nhets and age at collection - defunct
# plt.nhets.agec<-ggplot(hq.ndenovo.per.ind,aes(age.c,nhets))+
#   geom_point()+
#   geom_line(data=pred.age.c,aes(x=age.c,pred.nhets),color="red")+
#   theme_bw()+
#   labs(x="Age at collection of individual",y="No. of heteroplasmies in individual's tissues")+
#   annotate(geom="text",x=30,y=2.5,label="Beta= -0.010,\np-value= 0.611",size=5)

#plot relationship between nhets and age at birth
plt.nhets.ageb<-ggplot(hq.ndenovo.per.ind,aes(age.b,nhets))+
  geom_point()+
  geom_line(data=pred.age.b,aes(x=age.b,pred.nhets),color="red")+
  theme_bw()+
  labs(x="Age at collection of individual",y="No. of heteroplasmies in individual's tissues")+
  annotate(geom="text",x=25,y=2.5,label="Beta= 0.028,\np-value= 0.174",size=5)

ggsave("files/Analysis/Denovo_mutations/plt_ndenovo_hets_by_age_collection.pdf",plt.nhets.agec,height=5,width=5)
ggsave("files/Analysis/Denovo_mutations//plt_ndenovo_hets_by_age_birth.pdf",plt.nhets.ageb,height=5,width=5)

write.table(allfam.denovo.red,"files/Analysis/Denovo_mutations/allfam_denovo_03152019.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(pred.age.b,"files/Analysis/Denovo_mutations/hq_ndenovo_per_ind_03152019.txt",sep="\t",col.names=T,row.names=F,quote=F)

