#Author: AAZaidi

#code to explore pathogenicity of heteroplasmies
#Q. at what frequency do pathogenic variants exist in tissues

library(ggplot2)
library(dplyr)
library(data.table)
library(Biostrings)

#load heteroplasmies with sequence context info
hq.context<-fread("files/Analysis/syn_nsyn_annotation/hq_adjf_ancestral_syn_11282018.txt",header=T,sep="\t")

#make column for mutant frequency
#if the mutant frequency is the major allele, the mutant frequency is 1-MAF
hq.context2<-hq.context%>%
  filter(is.na(ancestral)=='FALSE'&ancestral.py!="")%>%
  mutate(mutant.f=case_when(major==ancestral~maf,
                            minor==ancestral~1-maf),
         mutation.id=paste(ancestral.py,position,mutant.py,sep=""))


#load pathogenicity annotation
#apogee+mitomap
pat.annot<-fread("files/Analysis/syn_nsyn_annotation/pathogenic_annotation_03222019.txt",header=T,sep="\t")

#add mutation_id column for merging with hq data
pat.annot<-pat.annot%>%
  mutate(mutation_id=paste(ancestral_py,position,mutant_py,sep=""))

#merge heteroplasmy data with pathogenicity annotation
hq.pat<-merge(hq.context2,pat.annot,by.x="mutation.id",by.y="mutation_id",all.x=T)

#add column for pathogenicity
hq.pat<-hq.pat%>%
  mutate(pathogenicity=case_when(is.na(source)=="FALSE" ~"P",
                                 is.na(source)=="TRUE" ~ "N"))

#calculate the median frequency for pathogenic vs non-pathpgenic variants
hqpat_freq<-hq.pat%>%
  group_by(tissue,pathogenicity)%>%
  summarize(median=median(mutant.f),mean=mean(mutant.f))

#test whether the mutant frequency is different between neutral and pathogenic variants
wtest.pat<-wilcox.test(data=hq.pat,mutant.f~pathogenicity)
wtest.pat.p<-formatC(wtest.pat$p.value,digits=2,format="e")

#test for difference in frequency b/w blood
wtest.pat.bl<-wilcox.test(data=hq.pat[which(hq.pat$tissue=="bl"),],mutant.f~pathogenicity)
wtest.pat.bl.p<-formatC(wtest.pat.bl$p.value,digits=2,format="e")

#test for difference in frequency b/w cheek
wtest.pat.ch<-wilcox.test(data=hq.pat[which(hq.pat$tissue=="ch"),],mutant.f~pathogenicity)
wtest.pat.ch.p<-formatC(wtest.pat.ch$p.value,digits=2,format="e")

#plot the distribution of mutant frequency for pat/non-pat variants
plt_hqpat<-ggplot()+
  geom_boxplot(data=hq.pat,aes(tissue,
                             mutant.f,
                             color=pathogenicity))+
  geom_point(data=hq.pat,aes(x=tissue,
                                 y=mutant.f,
                             color=pathogenicity),
             position=position_jitterdodge(jitter.width = 0.2),
             alpha=0.4)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size=14))+
  labs(x="Tissue type",
       y="Mutant allele frequency",
       color="Path.")+
  scale_fill_discrete(labels=c("N","P"))+
  scale_x_discrete(labels=c("Blood","Cheek"))+
  scale_y_log10(limits=c(0.005,1),
                breaks=c(0.01,seq(0.1,1,0.2)))

ggsave("files/Analysis/mutation_spectrum/plt_hqpat_f_04022019.pdf",plt_hqpat,height=7,width=7,useDingbats=F)



fwrite(hq.pat,"files/Analysis/syn_nsyn_annotation//hq_adjf_ancestral_syn_pathogenic_03222019.txt",col.names = T,row.names = F,quote=F,sep="\t")


#frequency b/w synonymous and non-synonymous mutations
hq.syn<-hq.context2%>%
  filter(syn!="")
  
hq.syn.summary<-hq.syn%>%
  group_by(tissue,syn)%>%
  summarize(median_f=median(mutant.f))

#test whether the mutant frequency is different between neutral and pathogenic variants
wtest.syn<-wilcox.test(data=hq.syn,mutant.f~syn)
wtest.syn.p<-formatC(wtest.syn$p.value,digits=2,format="e")


#test whether the mutant frequency is different in blood
wtest.syn.bl<-wilcox.test(data=hq.syn[which(hq.syn$tissue=="bl"),],mutant.f~syn)
wtest.syn.bl.p<-formatC(wtest.syn.bl$p.value,digits=2,format="e")

#test whether the mutant frequency is different in cheek
wtest.syn.ch<-wilcox.test(data=hq.syn[which(hq.syn$tissue=="ch"),],mutant.f~syn)
wtest.syn.ch.p<-formatC(wtest.syn.ch$p.value,digits=2,format="e")

hq.syn$syn<-factor(hq.syn$syn,levels=c("syn","nsyn"))

#plot the distribution of mutant frequency for syn/non-syn variants
plt_hqsyn<-ggplot()+
  geom_boxplot(data=hq.syn,aes(syn,
                             mutant.f,
                             color=tissue))+
  geom_point(data=hq.syn,aes(x=syn,
                             y=mutant.f,
                             color=tissue),
             position=position_jitterdodge(jitter.width=0.2),
             alpha=0.4)+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size=14))+
  labs(x="Tissue type",
       y="Mutant allele frequency",
       color="Syn/\nnonsyn")+
  scale_color_discrete(labels=c("S","N"))+
  scale_x_discrete(labels=c("Blood","Cheek"))+
  scale_y_log10(limits=c(0.005,1),
                breaks=c(0.01,seq(0.1,1,0.2)))

ggsave("files/Analysis/mutation_spectrum/plt_hqsyn_f_04022019.pdf",plt_hqsyn,
       height=7,
       width=7,useDingbats=F)


