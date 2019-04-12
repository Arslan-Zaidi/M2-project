#Author: AAZaidi

library(ggplot2)
library(dplyr)
library(data.table)
library(pbapply)
library(ape)
library(Biostrings)

#load rcrs fasta file
rcrs<-read.dna("../../../files/Database/rcrs.fasta",format="fasta")

#load mt genetic code
mt.code<-getGeneticCode("2")

#load genbank file - rcrs annotation for genes
genedb<-fread("../../../files/Database/genedb.txt",header=T)

#keep coding sequences only
genedb2<-genedb%>%
  filter(class=="CDS")

#calculate no. of syn and nysn sites for Nei-Gojobori method
genedb2$syn.sites<-NA
genedb2$nsyn.sites<-NA


#write function to calculate nsyn sites and syn sites for a single codon
nei.sites<-function(codon){
  aa<-translate(DNAString(codon),mt.code)
  codon<-unlist(strsplit(codon,split=""))
  s<-1
  for(k in 1:3){
    nucs<-setdiff(c("a","c","g","t"),codon[k])
    for(l in 1:3){
      mut.codon<-codon
      mut.codon[k]<-nucs[l]
      mt.aa<-translate(DNAString(paste(mut.codon,collapse="")),mt.code)
      if(mt.aa==aa){s<-s+1/3}else{s<-s}
    }
  }
  n<-3-s
  return(c(syn.sites=s,nsyn.sites=n))
}

#write function to return fold degeneracy for a site in a given codon
fold.degeneracy<-function(codon,site){
  aa<-translate(DNAString(codon),mt.code)
  codon<-unlist(strsplit(codon,split=""))
  nucs<-setdiff(c("a","c","g","t"),codon[site])
  s<-1
  for(l in 1:3){
    mut.codon<-codon
    mut.codon[site]<-nucs[l]
    mt.aa<-translate(DNAString(paste(mut.codon,collapse="")),mt.code)
    if(mt.aa==aa){s<-s+1}else{s<-s}
  }
  #if all 3 substitutions are nonsynonymous, site is 4-fold denegerate
  return(s)
}

#apply this function over each codon of every gene
genedb2$syn.sites<-NA
genedb2$nsyn.sites<-NA
genedb2$fold1<-NA
genedb2$fold2<-NA
genedb2$fold3<-NA
genedb2$fold4<-NA
for(i in 1:nrow(genedb2)){
  gene.start<-genedb2$start[i]+1
  gene.end<-genedb2$end[i]
  gene.seq<-as.character(rcrs[c(gene.start:gene.end)])
  if(genedb2$strand[i]==1){gene.seq<-gene.seq}
  if(genedb2$strand[i]==-1){
    gene.seq<-tolower(unlist(strsplit(as.character(reverseComplement(DNAString(paste(gene.seq,collapse="")))),split="")))}
  if(length(gene.seq)%%3==0){
    gene.seq<-gene.seq
  }
  if(length(gene.seq)%%3==1){
    gene.seq<-c(gene.seq,"a")
  }
  if(length(gene.seq)%%3==2){
    gene.seq<-c(gene.seq,"a","a")
  }
  gene.codons<-seqinr::splitseq(gene.seq,0,3)
  print(genedb2$gene[i])
  fold1<-0
  fold2<-0
  fold3<-0
  fold4<-0
  for(j in 1:length(gene.codons)){
    for(k in 1:3){
      f.d<-fold.degeneracy(gene.codons[j],k)
      if(f.d==1){fold1<-fold1+1}
      if(f.d==2){fold2<-fold2+1}
      if(f.d==3){fold3<-fold3+1}
      if(f.d==4){fold4<-fold4+1}
    }
  }
  genedb2$fold1[i]<-fold1
  genedb2$fold2[i]<-fold2
  genedb2$fold3[i]<-fold3
  genedb2$fold4[i]<-fold4
  
  nsites<-pbsapply(gene.codons,nei.sites)
  nsites<-apply(nsites,1,sum)
  genedb2$syn.sites[i]<-nsites[1]
  genedb2$nsyn.sites[i]<-nsites[2]
}

#syn/nsyn annotation
hq.syn<-fread("../../../files/Analysis/syn_nsyn_annotation/hq_syn_nsyn_annotation_10152018.txt",header=T,sep="\t")

#context annotation
hq.context<-fread("../../../files/Analysis/syn_nsyn_annotation/hq_conservative_ancestral_annotation_10212018.txt",header = T,sep="\t")

#how many independent mutations per mutation type
hq.context%>%
  filter(duplicated(fam_het_id)=="FALSE")%>%
  group_by(mutation)%>%
  summarize(n=length(maf))

#merge both
hq.merged<-hq.context%>%
  select(fqid,position,FID,mother_id,individual_id,level,tissue,tissue_id,Sex,haplogroup,age_collection,age_birth,fq_het_id,fam_het_id,tissue_het_id,ancestral,mutant,context,ancestral.py,mutant.py,context.py,mutation)%>%
  merge(.,hq.syn,by=c("fqid","position"))

#plot distribution of mutations along mtDNA
hq.fid.uniq<-hq.merged%>%
  filter(duplicated(fq_het_id)=="FALSE")%>%
  group_by(position)%>%
  summarize(nhets=length(unique(FID)))

plt.distr.hets<-ggplot()+
  geom_point(data=hq.fid.uniq,aes(position,nhets),size=0.5)+
  xlim(c(1,16569))+
  theme_bw()

##### plot proportion of heteroplasmies per fold degeneracy
dhq.fold<-hq.merged%>%
  filter(duplicated(fq_het_id)=="FALSE" & gene%in%genedb2$gene)%>%
  group_by(fold.d)%>%
  summarize(nhets=length(unique(FID)))

dgenedb.fold<-genedb2%>%
  select(gene,fold1,fold2,fold3,fold4)%>%
  melt(id.var="gene",variable.name="degeneracy",value.name="fold.d")%>%
  group_by(degeneracy)%>%
  summarize(sum.fold.d=sum(fold.d))

dmerged.fold<-cbind(hq.fold,genedb.fold)
dmerged.fold%>%
  mutate(prop.hets=nhets/sum.fold.d)%>%
  ggplot(.)+
  geom_bar(aes(fold.d,prop.hets),stat="identity")+
  theme_bw()+
  labs(x="Fold degeneracy",y="No. of hets/No. of sites")


###### plot het data per syn/nsyn
dhq.syn<-hq.merged%>%
  filter(duplicated(fq_het_id)=="FALSE" & gene%in%genedb2$gene)%>%
  group_by(gene,syn)%>%
  summarize(nhets=length(FID))

dgenedb.syn<-genedb2%>%
  mutate(syn=syn.sites,nsyn=nsyn.sites)%>%
  select(gene,syn,nsyn)%>%
  melt(id.var="gene",variable.name="synonymity",value.name="nsites")%>%
  group_by(gene,synonymity)%>%
  summarize(sum.nsites=sum(nsites))

dmerged.syn<-merge(dhq.syn,dgenedb.syn,by.x=c("gene","syn"),by.y=c("gene","synonymity"))

dmerged.syn%>%
  mutate(prop.hets=nhets/sum.nsites)%>%
  ggplot(.)+
  geom_bar(aes(syn,prop.hets),stat="identity")+
  facet_wrap(~gene)


#how many heteroplasmies for males and females
hq%>%
  group_by(Sex)%>%
  summarize(n=length(unique(individual_id)),nhets=length(unique(c(individual_id,position))))%>%
  mutate(prop=nhets/n)%>%
  ggplot()+
  geom_bar(aes(Sex,prop),stat="identity")


