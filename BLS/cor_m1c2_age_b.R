#Author: AAZaidi

library(ggplot2)
library(dplyr)
library(data.table)
library(purrr)

#### load heteroplasmies - conservative (one mother and two children families only) #####
hq.counts<-fread("files/Database/Heteroplasmy_tables/hq_counts_adj_frequency_09272018.txt",
                 header=T)

#add unique identifiers to table
hq.counts<-hq.counts%>%
  mutate(individual_het_id=paste(hq.counts$individual_id,hq.counts$position,sep="_"),
         mother_het_id=paste(hq.counts$mother_id,hq.counts$position,sep="_"))

#select specific columns to keep. We don't need all of them.
hq<-hq.counts%>%
  select(FID,fam_cat,fam_het_id,mother_id,mot_cat,mother_het_id,individual_id,individual_het_id,tissue_id,position,level,tissue,major,minor,adj.f,cvrg,age_collection,age_birth)


#divide data into two based on age at collection
q.age_collection<-hq%>%
  select(individual_id,age_collection)%>%
  distinct()%>%
  pull(age_collection)%>%
  quantile(.,probs=c(0,0.5,1))/365

#median age = 11.85

#mothers
hq.mothers<-hq%>%
  filter(mot_cat!="")%>%
  mutate(mother_group=individual_het_id,
         age_birth_cat=NA,
         age_collect_cat=cut(age_collection/365,breaks=q.age_collection,ordered_result=T,include.lowest=T))

#their children
hq.children<-hq%>%
  filter(mother_id%in%hq.mothers$individual_id)

#get age of mother at childbirth for these children
q.age_birth<-hq.children%>%
  select(mother_id,age_birth)%>%
  distinct()%>%
  pull(age_birth)%>%
  quantile(.,probs=c(0,0.5,1))/365

hq.children<-hq.children%>%
  mutate(mother_group=mother_het_id,
         age_birth_cat=cut(age_birth/365,breaks=q.age_birth,ordered_result=T,include.lowest=T),
         age_collect_cat=cut(age_collection/365,q.age_collection,ordered_result=T,include.lowest=T))

#average frequency across bl and ch tissues
hq.mothers<-hq.mothers%>%
  select(individual_id,position,individual_het_id,tissue,adj.f)%>%
  group_by(individual_id,position,individual_het_id)%>%
  summarize(adj.f=mean(adj.f))

#calculate average frequency b/w bl and ch tissues of children
hq.children<-hq.children%>%
  select(mother_id,individual_id,position,mother_het_id,individual_het_id,tissue,adj.f,age_birth_cat)%>%
  group_by(mother_id,mother_het_id,individual_id,position,individual_het_id,age_birth_cat)%>%
  summarize(adj.f=mean(adj.f))%>%
  select(mother_id,position,mother_het_id,adj.f,age_birth_cat)

#merge mother and child data
hq.merged<-merge(hq.mothers,hq.children,by.x="individual_het_id",by.y="mother_het_id")

#calculate correlation b/w allele frequency between mothers and children
cor.summary<-hq.merged%>%
  group_by(age_birth_cat)%>%
  summarize(correlation=cor(adj.f.x,adj.f.y))%>%
  mutate(y=c(0.9,1))

mc_bage_cor<-ggplot()+
  geom_point(data=hq.merged,aes(adj.f.x,adj.f.y,color=age_birth_cat),alpha=0.5)+
  theme_bw()+
  scale_color_manual(values=c("#1b9e77","#d95f02"),
                     labels=c("[15-29]","[29-46]"))+
  geom_text(data=cor.summary,aes(x=0.1,y=y,label=paste("r=",round(correlation,2)),color=age_birth_cat))+
  labs(x="Heteroplasmy allele frequency - mother",
       y="Heteroplasmy allele frequency - child",
       color="Age of mother at\nchild birth")

ggsave("files/Analysis/Branch_stats/plt_m1c_ageb_cor_03202019.pdf",width=7,height=5,useDingbats=F)
