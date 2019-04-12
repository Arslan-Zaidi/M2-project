library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)

#load genbank file - rcrs annotation for genes
genedb<-fread("files/Database/genedb.txt",header=T)
genedb<-genedb%>%
  mutate(region=case_when(startsWith(gene,"TRN")~"TRNA",
                   startsWith(gene,"RNR")~"RNA",
                   startsWith(gene,"D-loop")~"D-loop",
                   TRUE~"Coding"))

genedb$gene<-factor(genedb$gene,levels=genedb$gene[order(-genedb$start2)])


hq.context<-fread("files/Analysis/syn_nsyn_annotation/hq_adjf_ancestral_syn_pathogenic_11282018.txt",header=T,sep="\t")
hq2<-hq.context%>%
  distinct(fam_het_id,.keep_all=T)

hq2<-hq2%>%
  mutate(mutation=paste(major,minor,sep=""))%>%
  mutate(tstv=ifelse(mutation%in%c("AG","GA","CT","TC"),"ts","tv"),
         syn=ifelse(syn=="","non-coding",syn),
         region=case_when(startsWith(gene,"TRN")~"TRNA",
                          startsWith(gene,"RNR")~"RNA",
                          startsWith(gene,"D-loop")~"D-loop",
                          TRUE~"Coding"))

p1<-ggplot(hq2)+
  geom_point(alpha=0.8,size=1.5,aes(position,FID,shape=tstv,color=region))+
  theme_bw()+
  theme(axis.text.y=element_blank(),panel.grid.major.y = element_blank())+
  scale_shape_manual(values=c(1,19))+
  scale_color_manual(values=c("#984ea3","#e41a1c","#377eb8","#ff7f00"))+
  labs(y="Families",x="rCRS coordinates",color="Region",shape="Ts/Tv")
p1

ggsave("files/Analysis/mutation_spectrum/plt_hq_alongmtdna_03182019.pdf",p1,height=5,width=7)

dhq<-hq2%>%
  group_by(position)%>%
  summarize(nhets=length(unique(FID)))

p2<-ggplot(dhq,aes(position,nhets))+
  geom_line()+
  theme_void()
p2

dhq3<-hq2%>%
  group_by(FID)%>%
  summarise(nhets=length(unique(position)))

p3<-ggplot(dhq3,aes(FID,nhets))+
  geom_bar(stat="identity",width=1)+
  theme_classic()+theme(axis.text.x=element_blank(),
                        axis.title.x = element_blank(),
                        axis.line.x = element_blank())+
  scale_y_continuous(breaks=seq(0,14,2))+
  coord_flip()

p3

p2<-ggplot(genedb)+
  geom_segment(aes(x=start2,xend=end2,y=gene,yend=gene,color=region),arrow=arrow(length=unit(0.05,"inches")))+
  theme_bw()+
  scale_color_manual(values=c("#984ea3","#e41a1c","#377eb8","#ff7f00"))+
  labs(x=)
p2  

