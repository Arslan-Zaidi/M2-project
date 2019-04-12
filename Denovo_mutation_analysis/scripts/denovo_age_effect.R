
library(dplyr)
library(ggplot2)

allfam.denovo<-read.table("../../../files/Analysis/Denovo_mutations/allfam_denovo_10122018.txt",header=T,sep="\t")

hq.denovo<-allfam.denovo%>%
  filter(maf>0.01)

ddply(hq.denovo,.(FID,individual_id,position),summarize,n=length(maf))

#explore the effects of age on no. of denovo heteroplasmies
hq.denovo.age<-hq.denovo%>%
  group_by(individual_id,age_collection,age_birth)%>%
  summarize(nhets=length(unique(position)))

hq.denovo.age<-hq.denovo.age%>%
  mutate(age.b=age_birth/365,age.c=age_collection/365)

#poisson regression on age at collection
hq.age.c<-glm(data=hq.denovo.age,nhets~age.c,family="poisson")
hq.denovo.age$nhets.c<-predict(hq.age.c,type="response")

#poisson regression on age at birth
hq.age.b<-glm(data=hq.denovo.age,nhets~age.b,family="poisson")
hq.denovo.age$nhets.b<-predict(hq.age.b,type="response")

plt.nhets_agec<-ggplot(hq.denovo.age,aes(age.c,nhets))+
  geom_point()+
  geom_line(aes(x=age.c,nhets.c),color="red")+
  theme_bw()+
  labs(x="Age at collection of individual",y="No. of heteroplasmies in individual's tissues")+
  annotate(geom="text",x=25,y=3,label="Beta= -0.01,\np-value= 0.46",size=5)

plt.nhets_ageb<-ggplot(hq.denovo.age,aes(age.b,nhets))+geom_point()+
  geom_line(aes(x=age.b,y=nhets.b),color="red")+
  theme_bw()+
  labs(x="Age of mother at birth of child",y="No. of heteroplasmies in child")+
  annotate(geom="text",x=25,y=3,label="Beta= 0.1,\np-value= 0.74",size=5)

ggsave("../../../files/Analysis/PNAS_styled_analyses/plt_ndenovo_hets_by_age_collection.pdf",plt.nhets_agec,height=5,width=5)
ggsave("../../../files/Analysis/PNAS_styled_analyses/plt_ndenovo_hets_by_age_birth.pdf",plt.nhets_ageb,height=5,width=5)


# test for non-uniform distribution of de novo mutations

hq.denovo<-merge(hq.denovo,annot,by="fq_het_id")
dhq.denovo<-hq.denovo%>%
  group_by(gene)%>%
  summarize(nhets=length(unique(individual_id)))

dhq.denovo<-merge(dhq.denovo,genedb[c('start','end','class','gene')],by="gene")  
dhq.denovo$length<-with(dhq.denovo,end-start)

#expected number of hets per gene
dhq.denovo$exp.nhets<-NA
dhq.denovo$exp.nhets<-with(dhq.denovo,sum(nhets)/sum(length) * length)
ddhq.denovo<-dhq.denovo%>%
  group_by(class)%>%
  summarize(obs.nhets=sum(nhets),exp.nhets=sum(exp.nhets))
mddhq.denovo<-melt(ddhq.denovo,id.vars="class")

#compare this with distribution of ALL heteroplasmies
hq.all<-fread("../../../files/Analysis/Fixing_heteroplasmy_table/hq_cleaned_conservative_09272018.txt",header=T,sep="\t")

hq.all<-merge(hq.all,annot,by="fq_het_id")
dhq.all<-hq.all%>%
  group_by(gene)%>%
  summarize(nhets=length(unique(fam_het_id)))


dhq.all<-merge(dhq.all,genedb[c('start','end','class','gene')],by="gene")  
dhq.all$length<-with(dhq.all,end-start)

#expected number of hets per gene
dhq.all$exp.nhets<-NA
dhq.all$exp.nhets<-with(dhq.all,sum(nhets)/sum(length) * length)
ddhq.all<-dhq.all%>%
  group_by(class)%>%
  summarize(obs.nhets=sum(nhets),exp.nhets=sum(exp.nhets))

mddhq.all<-melt(ddhq.all,id.vars="class")

mddhq.all$heteroplasmy_type<-"All"
mddhq.denovo$heteroplasmy_type="De Novo"

mddhq.combined<-rbind(mddhq.all,mddhq.denovo)

plt.distribution.by.class<-ggplot(mddhq.combined,aes(x=class,y=value))+
  geom_bar(aes(fill=variable),stat="identity",position="dodge")+
  facet_wrap(~heteroplasmy_type,scale="free")+
  labs(x="Functional class",y="Counts",fill="Expected/observed\nno. of heteroplasmies")+
  theme_bw()

ggsave("../../../files/Analysis/PNAS_styled_analyses/plt_hets_distr_by_class.pdf",plt.distribution.by.class,height=5,width=10)

#apply chi.sq test to regions
#write function
chi.test<-function(x){
  #x is the table of obs and expected values
  
}



#explore the effects of age on no. of denovo heteroplasmies
hq.all.age<-hq.all%>%
  group_by(individual_id,age_collection,age_birth)%>%
  summarize(nhets=length(unique(position)))

hq.all.age<-hq.all.age%>%
  mutate(age.b=age_birth/365,age.c=age_collection/365)

#poisson regression on age at collection
hq.age.c<-glm(data=hq.all.age,nhets~age.c,family="poisson")
hq.all.age$nhets.c<-predict(hq.age.c,type="response")

#poisson regression on age at birth
hq.age.red<-hq.all.age%>%filter(is.na(age.b)=='FALSE')
hq.age.b<-glm(data=hq.age.red,nhets~age.b,family="poisson")
hq.age.red$nhets.b<-predict(hq.age.b,type="response")

plt.nhets_agec<-ggplot(hq.all.age,aes(age.c,nhets))+
  geom_point()+
  geom_line(aes(x=age.c,nhets.c),color="red")+
  theme_bw()+
  labs(x="Age at collection of individual",y="No. of heteroplasmies in individual's tissues")+
  annotate(geom="text",x=50,y=7,label="Beta= 0.002,\np-value= 0.45",size=5)

plt.nhets_ageb<-ggplot(hq.age.red,aes(age.b,nhets))+geom_point()+
  geom_line(aes(x=age.b,y=nhets.b),color="red")+
  theme_bw()+
  labs(x="Age of mother at birth of child",y="No. of heteroplasmies in child")+
  annotate(geom="text",x=30,y=7,label="Beta= 012,\np-value= 0.14",size=5)

ggsave("../../../files/Analysis/PNAS_styled_analyses/plt_all_hets_by_age_collection.pdf",plt.nhets_agec,height=5,width=5)
ggsave("../../../files/Analysis/PNAS_styled_analyses/plt_all_hets_by_age_birth.pdf",plt.nhets_ageb,height=5,width=5)



