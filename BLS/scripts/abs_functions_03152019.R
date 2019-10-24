#Author: AAZaidi

########### functions to calculate Fst components and branch lengths ######################
library(dplyr)

#function to calculate the numerator (num) and denominator (den) for fst calculation
#fst formulae from Bahatia et al. 2013 Gen. Research and Reynold's 1983
fst_comps<-function(p1,p2,n1,n2,estimator="h"){
  #p1 is allele frequency in pop1
  #p2 is allele frequency in pop2
  #n1 is sample size of pop1
  #n2 is sample size of pop2
  if(!estimator%in%c("r","h","wc","b")){stop("estimator can either be 'h' (hudson),'r' (reynolds), 'wc'(weir-cockerham), or 'b' (basic)")}
  if(estimator=="r"){
    h1=2*p1*(1-p1)
    h2=2*p2*(1-p2)
    num=(p1-p2)^2 - (n1+n2)*(n1*h1 + n2*h2)/(4*n1*n2*(n1+n2 -1))
    den=(p1-p2)^2 + ((4*n1*n2 - n1 - n2)*(n1*h1+n2*h2)/(4*n1*n2*(n1+n2 -1)))
  }
  if(estimator=="h"){
    num=(p1-p2)^2 - ( (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1) )
    den=(p1*(1-p2)) + (p2*(1-p1))
  }
  if(estimator=="wc"){
    a = 2*(n1*n2)/(n1+n2) * 1/(n1+n2 -2) * (n1*(p1*(1-p1)) + n2*(p2*(1-p2)))
    b = (n1*n2)/(n1+n2) * (p1-p2)^2 + ((2*n1*n2/(n1+n2)) -1) * (1/(n1 + n2 -2)) * (n1*p1*(1-p1) + n2*(p2*(1-p2)))
    num=b-a
    den=b
  }
  if(estimator=="b"){
    a=(2*p1*(1-p1) + 2*p2*(1-p2))/2
    pbar=(p1+p2)/2
    b=2*pbar*(1-pbar)
    num=b-a
    den=b
  }
  return(c(num=num,den=den))
}


#function to calculation divergence between x and y (Dxy) from fst
dxy<-function(fst){
  if(is.na(fst)=="FALSE"){
  if(fst==1){
    d<- -2*log(1-0.99)}else{d<- -2*log(1-fst)}
    }else{d<-NA}
  return(d)}

cal_fst<-function(dat,est="h"){
  dat$lt.comb<-with(dat,paste(level,tissue,sep="_"))
  leaves<-sort(unique(dat$lt.comb))
  l.combs<-combn(leaves,2)
  l1<-l.combs[1,]
  l2<-l.combs[2,]
  fst.out<-mapply(function(x,y){
    dat1<-filter(dat,lt.comb==x)
    dat2<-filter(dat,lt.comb==y)
    info1<-dat1%>%select(.,adj.f,cvrg)
    info2<-dat2%>%select(.,adj.f,cvrg)
    
    freq1<-as.numeric(info1[1])
    cvrg1<-as.numeric(info1[2])
    
    freq2<-as.numeric(info2[1])
    cvrg2<-as.numeric(info2[2])
    output<-fst_comps(freq1,freq2,cvrg1,cvrg2,est)
    output<-as.data.frame(output)
    output$leaf1<-x
    output$leaf2<-y
    
    spltx1<-unlist(strsplit(x,split="_"))
    spltx2<-unlist(strsplit(y,split="_"))
    
    output$ind1<-spltx1[1]
    output$ind2<-spltx2[1]
    
    output$tis1<-spltx1[2]
    output$tis2<-spltx2[2]
    
    output$age.c1<-dat1%>%pull(age_collection)/365
    output$age.c2<-dat2%>%pull(age_collection)/365
    
    output$age.b1<-dat1%>%pull(age_birth)/365
    output$age.b2<-dat2%>%pull(age_birth)/365
    
    output$age.b1.cat<-dat1%>%pull(age_birth_cat)
    output$age.b2.cat<-dat2%>%pull(age_birth_cat)
    
    output$age.c1.cat<-dat1%>%pull(age_collect_cat)
    output$age.c2.cat<-dat2%>%pull(age_collect_cat)
    
    output$component<-row.names(output)
    return(output)
  },l1,l2,SIMPLIFY = F)
  names(fst.out)<-mapply(paste,sep=".",l1,l2)
  fst.out<-bind_rows(fst.out,.id="pair")
  #add family category information
  fst.out$fam_cat<-unique(dat$fam_cat)
    #add mother category information
    fst.out$mot_cat<-unique(
      dat$mot_cat[which(dat$mother_group==dat$individual_het_id)])
  #add mean.f across leaves of the tree
  #for filtering later on allele frequency
  #fst.out$max.f<-max(dat$adj.f)

  #number of heteroplasmic individuals (maf>0.01)
  fst.out$max.f<-length(which(dat$adj.f>0.01 & dat$adj.f<1 | 
                                dat$adj.f<0.99 & dat$adj.f>0))
  
  #add other annotations
  fst.out$region<-unique(dat$region)
  fst.out$gene<-unique(dat$gene)
  fst.out$fold.d<-unique(dat$fold.d)
  fst.out$syn<-unique(dat$syn)
  fst.out$pathogenicity<-unique(dat$pathogenicity)
  fst.out$mother_group=unique(dat$mother_group)
  return(fst.out)
}

#function to calculate BLS across all heteroplasmies
cal_bls<-function(fstlist=shq.fst,min.hets,mot_cats){
  if(mot_cats=="all"){
    fstlist2=fstlist%>%
      keep(function(x){mean(x$max.f)>=min.hets})
  }else{
  fstlist2=fstlist%>%
    keep(function(x){mean(x$max.f)>=min.hets & x$mot_cat[1]%in%mot_cats})
  }

  fstlist3<-bind_rows(unname(fstlist2),.id="sno")%>%
    mutate(branch.type=case_when(pair%in%c("c1_bl.c1_ch","c1_ch.c1_bl","m1_bl.m1_ch","m1_ch.m1_bl","c2_bl.c2_ch","c2_ch.c2_bl")~"bl2ch.div",
                                 pair%in%c("c1_bl.m1_bl","c2_bl.m1_bl")~"m_bl.c_bl",
                                 pair%in%c("c1_bl.m1_ch","c2_bl.m1_ch","c1_ch.m1_bl","c2_ch.m1_bl")~"m_bl.c_ch",
                                 pair%in%c("c1_ch.m1_ch","c2_ch.m1_ch")~"m_ch.c_ch",
                                 pair%in%c("c1_bl.c2_bl")~"c_bl.c_bl",
                                 pair%in%c("c1_ch.c2_ch")~"c_ch.c_ch",
                                 pair%in%c("c1_bl.c2_ch","c1_ch.c2_bl")~"c_bl.c_ch"))
  
    
  
  
  fstlist4<-fstlist3%>%
    dcast(sno+mother_group+pair+leaf1+leaf2+ind1+ind2+tis1+tis2+fam_cat+branch.type~component,value.var="output")
  
  #calculate fst
  #fst1: ratio of avgs
  fstlist5<-fstlist4%>%
    group_by(pair)%>%
    summarize(fst1=sum(num)/sum(den))
  
  #calculate dxy (-2*log(1-fst)) - this is the branch length
  fstlist5$dxy1<-sapply(fstlist5$fst1,dxy)
  fst.df<-fstlist5
  fst.df$sno<-1
  bls.df<-fst.df%>%
    dcast(sno~pair,value.var="dxy1")%>%
    mutate(m1_c1=0.5*max((c1_bl.m1_bl+c1_ch.m1_ch),(c1_bl.m1_ch+c1_ch.m1_bl)) - (m1_bl.m1_ch+c1_bl.c1_ch),
           m1_c2=0.5*max((c2_bl.m1_bl+c2_ch.m1_ch),(c2_bl.m1_ch+c2_ch.m1_bl)) - (m1_bl.m1_ch+c2_bl.c2_ch),
           c1_c2=0.5*max((c1_bl.c2_bl+c1_ch.c2_ch),(c1_bl.c2_ch+c1_ch.c2_bl)) - (c1_bl.c1_ch+c2_bl.c2_ch),
           m1_c=0.5*((m1_c1+m1_c2) - c1_c2),
           m.tdiv=m1_bl.m1_ch,
           c1.tdiv=c1_bl.c1_ch,
           c2.tdiv=c2_bl.c2_ch)
  bls.df$min.hets=min.hets
  bls.df$nfams=length(fstlist2)
  return(bls.df)
  
}

#function to calculate BLS for different age categories
#here we don't have to calculate internal branches
cal_bls_age_c<-function(fstlist,bootstrap=TRUE){
  if(bootstrap==TRUE){#bootstrap within each age_c category and calculate bls
    #convert list output to data.frame
    fstlist2<-bind_rows(unname(fstlist),.id="sno")
    #dcast so numerator and denominator are in different columns
    fstlist3<-fstlist2%>%
      filter(ind1==ind2)%>%
      dcast(sno+mother_group+pair+age.c1.cat~component,value.var="output")
    #calculate fst
    #fst1: ratio of avgs
    fstlist4<-fstlist3%>%
      group_by(age.c1.cat)%>%
      summarize(fst1=sum(num)/sum(den))
    
    #calculate dxy (-2*log(1-fst))
    fstlist4$dxy1<-sapply(fstlist4$fst1,dxy)
    
    #melt df to long format
    mfst.df<-fstlist4%>%
      select(age.c1.cat,dxy1)%>%
      rename(age_collection=age.c1.cat,length=dxy1)
  }
  if(bootstrap==FALSE){#average bls across heteroplasmies for each mother-child comparison
    #no bootstrapping
    #convert list output to data.frame
    fstlist2<-bind_rows(unname(fstlist),.id="sno")
    #dcast so numerator and denominator are in different columns
    fstlist3<-fstlist2%>%
      filter(ind1==ind2)%>%
      dcast(sno+mother_group+pair+age.c1~component,value.var="output")
    #calculate fst
    #fst1: ratio of avgs
    fstlist4<-fstlist3%>%
      group_by(age.c1)%>%
      summarize(fst1=sum(num)/sum(den))
    
    #calculate dxy (-2*log(1-fst))
    fstlist4$dxy1<-sapply(fstlist4$fst1,dxy)
    
    #melt df to long format
    mfst.df<-fstlist4%>%
      select(age.c1,dxy1)%>%
      rename(age_collection=age.c1,length=dxy1)
  }
  
  
  return(mfst.df)
  
}

#there are two functions for effect of age on mother's age at child birth
#1. bootstrap for each age_birth category and calculate distribution of bls
cal_bls_age_b1<-function(fstlist){
  #convert list output to data.frame
  fstlist2<-bind_rows(unname(fstlist),.id="sno")
  
  #dcast so numerator and denominator are in different columns
  fstlist3<-fstlist2%>%
    dcast(sno+mother_group+pair+leaf1+leaf2+ind1+ind2+tis1+tis2+fam_cat~component,value.var="output")
  
  #calculate fst
  #fst1: ratio of avgs
  fstlist4<-fstlist3%>%
    group_by(pair)%>%
    summarize(fst1=sum(num)/sum(den))
  
  #calculate dxy (-log(1-fst))
  fstlist4$dxy1<-sapply(fstlist4$fst1,dxy)

  
  fstlist4$sno<-1
  
  #calculate BLS using pairwise divergences
  bls.df<-fstlist4%>%
    dcast(sno~pair,value.var="dxy1")%>%
    mutate(m1_c1=0.5*max((c1_bl.m1_bl+c1_ch.m1_ch),(c1_bl.m1_ch+c1_ch.m1_bl)) - (m1_bl.m1_ch+c1_bl.c1_ch),
           m1_c2=0.5*max((c2_bl.m1_bl+c2_ch.m1_ch),(c2_bl.m1_ch+c2_ch.m1_bl)) - (m1_bl.m1_ch+c2_bl.c2_ch))
  
  #melt to long format for plotting etc.
  mbls.df<-bls.df%>%
    mutate(id="age_b")%>%
    select(id,m1_c1,m1_c2)%>%
    melt(id.var="id",variable.name="m_c",value.name="length")
  
  return(mbls.df)
}

#2. don't bootstrap. rather, calculate bls for each mother-child comparison, averaging across heteroplasmies
cal_bls_age_b2<-function(fstlist){
  #relabel mother-child pairs for consistency when averaging fst
  fstlist3<-lapply(fstlist,function(x){
    #mother category
    mot_category=unique(x$mot_cat)
    #dcast so numerator and denominator are in different columns
    fstlist2<-x%>%
      separate(mother_group,into=c("mother_id","position"),sep="_",remove=FALSE)%>%
      dcast(mother_id+position+pair+leaf1+leaf2+ind1+ind2+tis1+tis2+mot_cat+age.b1+age.b2~component,value.var="output")
    
    #relabel pair id so that there are only the following comparisons:
    ##m1_c1
    ##m1_m1
    ##c1_c1
    if(mot_category%in%c("m1c2","m1c3","m1c4","m1c5")){
      fstlist3<-fstlist2%>%
        mutate(pair2=case_when(ind1==ind2&ind1%in%c("c1","c2","c3","c4","c5")~paste("c1",tis1,"c1",tis2,sep="_"),
                               ind1==ind2&ind1%in%c("m1")~paste("m1",tis1,"m1",tis2,sep="_"),
                               ind1%in%c("c1","c2","c3","c4","c5")&ind2==("m1")~paste("c1",tis1,"m1",tis2,sep="_"),
                               
                               TRUE~"exclude"),
               age_birth=age.b1)
    }
    if(mot_category=="g1m2"){
      fstlist3<-fstlist2%>%
        mutate(pair2=case_when(ind1==ind2&ind1%in%c("g1")~paste("m1",tis1,"m1",tis2,sep="_"),
                               ind1==ind2&ind1%in%c("m1","m2")~paste("c1",tis1,"c1",tis2,sep="_"),
                               ind1%in%c("g1")&ind2%in%c("m1","m2")~paste("c1",tis2,"m1",tis1,sep="_"),
                               
                               TRUE~"exclude"),
               age_birth=age.b2)
    }
    if(mot_category=="t1g2"){
      fstlist3<-fstlist2%>%
        mutate(pair2=case_when(ind1==ind2&ind1%in%c("g1","g2")~paste("c1",tis1,"c1",tis2,sep="_"),
                               ind1==ind2&ind1%in%c("t")~paste("m1",tis1,"m1",tis2,sep="_"),
                               ind1%in%c("g1","g2")&ind2%in%c("t")~paste("c1",tis1,"m1",tis2,sep="_"),
                               
                               TRUE~"exclude"),
               age_birth=age.b1)
      
    }
    return(fstlist3)
  })
  fstlist4<-bind_rows(fstlist3)
  #calculate fst, divergence etc. for each mother-child pair
  fstlist4<-fstlist4%>%
    filter(pair2!="exclude")%>%
    group_by(mot_cat,mother_id,pair2,age_birth)%>%
    summarize(fst1=sum(num)/sum(den))
  #calculate dxy (-2log(1-fst))
  fstlist4$dxy1<-sapply(fstlist4$fst1,dxy)
  
  #split into list so that each mother-child transmission can be analyzed separately
  fstlist5<-split(fstlist4,fstlist4$mother_id)
  fstlist6<-lapply(fstlist5,function(x){
    #separately calculate divergence b/w mother's tissues
    #imp because mothers don't often have age of birth
    mtis_div<-x$dxy1[which(x$pair2=="m1_bl_m1_ch")]
    
    #calculate BLS using pairwise divergences
    bls.df<-x%>%
      filter(pair2!="m1_bl_m1_ch")%>%
      dcast(mother_id+age_birth~pair2,value.var="dxy1",fun.aggregate=mean)%>%
      mutate(m1_c1=0.5*max((c1_bl_m1_bl+c1_ch_m1_ch),(c1_bl_m1_ch+c1_ch_m1_bl)) - (mtis_div+c1_bl_c1_ch))%>%
      select(mother_id,age_birth,m1_c1)
    return(bls.df)
    
  })

  fstlist7<-bind_rows(fstlist6)
  
  return(fstlist7)
  
}

#function to calculate BLS for synonymous mutations
cal_bls_syn<-function(fstlist){
    
    #convert list output to data.frame
    boot.df<-bind_rows(unname(fstlist),.id="sno")
    
    boot.df2<-boot.df%>%
      dcast(sno+syn+mother_group+pair+leaf1+leaf2+ind1+ind2+tis1+tis2+fam_cat+branch.type~component,value.var="output")
    
    #calculate fst
    #fst2: avg. of ratios
    #fst1: ratio of avgs
    boot.df3<-boot.df2%>%
      group_by(syn,pair)%>%
      summarize(fst1=sum(num)/sum(den))
    
    #calculate dxy (-2*log(1-fst))
    boot.df3$dxy1<-sapply(boot.df3$fst1,dxy)
    
    #calculate BLS using pairwise divergences
    bls.boot<-boot.df3%>%
      dcast(syn~pair,value.var="dxy1")%>%
      mutate(m1_c1=0.5*(max((c1_bl.m1_bl+c1_ch.m1_ch),(c1_bl.m1_ch+c1_ch.m1_bl)) - (m1_bl.m1_ch+c1_bl.c1_ch)),
             m1_c2=0.5*(max((c2_bl.m1_bl+c2_ch.m1_ch),(c2_bl.m1_ch+c2_ch.m1_bl)) - (m1_bl.m1_ch+c2_bl.c2_ch)),
             c1_c2=0.5*(max((c1_bl.c2_bl+c1_ch.c2_ch),(c1_bl.c2_ch+c1_ch.c2_bl)) - (c1_bl.c1_ch+c2_bl.c2_ch)),
             m1_c=0.5*((m1_c1+m1_c2) - c1_c2),
             m.tdiv=m1_bl.m1_ch,
             c1.tdiv=c1_bl.c1_ch,
             c2.tdiv=c2_bl.c2_ch)
    
    #melt to long format for plotting etc.
    #melt the df to long format
    mbls.boot<-bls.boot%>%
      select(syn,m1_c1,m1_c2,c1_c2,m1_c,m.tdiv,c1.tdiv,c2.tdiv)%>%
      melt(id.var=c("syn"),value.name="length",variable.name="branch")%>%
      mutate(branch.type=case_when(branch=="m1_c1"~"m1_child1",
                                   branch=="m1_c2"~"m1_child2",
                                   branch=="c1_c2"~"c1_c2",
                                   branch=="m1_c"~"m1_children",
                                   branch=="m.tdiv"~"tissue.div",
                                   branch=="c1.tdiv"~"tissue.div",
                                   branch=="c2.tdiv"~"tissue.div"))
    
    mbls.boot$branch.type<-factor(mbls.boot$branch.type,
                                    levels=c("tissue.div","m1_child1","m1_child2","c1_c2","m1_children"))
    
    return(mbls.boot)
}

#calculate bls for pathogenic mutations
cal_bls_pat<-function(fstlist){
  
  #convert list output to data.frame
  boot.df<-bind_rows(unname(fstlist),.id="sno")
  
  boot.df2<-boot.df%>%
    dcast(sno+pathogenicity+mother_group+pair+leaf1+leaf2+ind1+ind2+tis1+tis2+fam_cat+branch.type~component,value.var="output")
  
  #calculate fst
  #fst2: avg. of ratios
  #fst1: ratio of avgs
  boot.df3<-boot.df2%>%
    group_by(pathogenicity,pair)%>%
    summarize(fst1=sum(num)/sum(den))
  
  #calculate dxy (-2*log(1-fst))
  boot.df3$dxy1<-sapply(boot.df3$fst1,dxy)
  
  #calculate BLS using pairwise divergences
  bls.boot<-boot.df3%>%
    dcast(pathogenicity~pair,value.var="dxy1")%>%
    mutate(m1_c1=0.5*(max((c1_bl.m1_bl+c1_ch.m1_ch),(c1_bl.m1_ch+c1_ch.m1_bl)) - (m1_bl.m1_ch+c1_bl.c1_ch)),
           m1_c2=0.5*(max((c2_bl.m1_bl+c2_ch.m1_ch),(c2_bl.m1_ch+c2_ch.m1_bl)) - (m1_bl.m1_ch+c2_bl.c2_ch)),
           c1_c2=0.5*(max((c1_bl.c2_bl+c1_ch.c2_ch),(c1_bl.c2_ch+c1_ch.c2_bl)) - (c1_bl.c1_ch+c2_bl.c2_ch)),
           m1_c=0.5*((m1_c1+m1_c2) - c1_c2),
           m.tdiv=m1_bl.m1_ch,
           c1.tdiv=c1_bl.c1_ch,
           c2.tdiv=c2_bl.c2_ch)
  
  #melt to long format for plotting etc.
  #melt the df to long format
  mbls.boot<-bls.boot%>%
    select(pathogenicity,m1_c1,m1_c2,c1_c2,m1_c,m.tdiv,c1.tdiv,c2.tdiv)%>%
    melt(id.var=c("pathogenicity"),value.name="length",variable.name="branch")%>%
    mutate(branch.type=case_when(branch=="m1_c1"~"m1_child1",
                                 branch=="m1_c2"~"m1_child2",
                                 branch=="c1_c2"~"c1_c2",
                                 branch=="m1_c"~"m1_children",
                                 branch=="m.tdiv"~"tissue.div",
                                 branch=="c1.tdiv"~"tissue.div",
                                 branch=="c2.tdiv"~"tissue.div"))
  
  mbls.boot$branch.type<-factor(mbls.boot$branch.type,
                                levels=c("tissue.div","m1_child1","m1_child2","c1_c2","m1_children"))
  
  return(mbls.boot)
  
}

#write function to sample 2 kids from each mother (if she has more than 2)
#also reformat dataframe for bls calculations
#this is useful if fst needs to be calculated using the ratio of averages method
sample2kids<-function(fst.obj){
  mot_cat=unique(fst.obj$mot_cat)
  #if 2 kids, return without changing
  if(mot_cat==c("m1c2")){fst.obj.tmp<-fst.obj}
  #if>2kids, sample 1 and rename them to be c1 and c2 based on age
  if(mot_cat%in%c("m1c3","m1c4","m1c5")){
    children.levels<-setdiff(fst.obj$ind1,"m1")
    sampled.children<-sort(sample(children.levels,2,replace=F))
    fst.obj.tmp<-fst.obj[which(fst.obj$ind1%in%c("m1",sampled.children) &
                                 fst.obj$ind2%in%c("m1",sampled.children)),]
    fst.obj.tmp<-fst.obj.tmp%>%
      mutate(ind1=case_when(ind1==sampled.children[1]~"c1",
                            ind1==sampled.children[2]~"c2",
                            !ind1%in%sampled.children~ind1),
             ind2=case_when(ind2==sampled.children[1]~"c1",
                            ind2==sampled.children[2]~"c2",
                            !ind2%in%sampled.children~ind2))%>%
      mutate(pair=paste(ind1,"_",tis1,".",ind2,"_",tis2,sep=""))
  }
  #standardize notation for g1m2 and t1g2
  if(mot_cat%in%c("g1m2")){
    fst.obj.tmp<-fst.obj%>%
      mutate(ind1=case_when(ind1=="m1"~"c1",
                            ind1=="m2"~"c2",
                            ind1=="g1"~"m1"),
             ind2=case_when(ind2=="m1"~"c1",
                            ind2=="m2"~"c2",
                            ind2=="g1"~"m1"))%>%
      mutate(pair=case_when(ind1=="c1"&ind2=="c2"~paste(ind1,"_",tis1,".",ind2,"_",tis2,sep=""),
                            ind1=="m1"&ind2=="c1"~paste(ind2,"_",tis2,".",ind1,"_",tis1,sep=""),
                            ind1=="m1"&ind2=="c2"~paste(ind2,"_",tis2,".",ind1,"_",tis1,sep=""),
                            ind1==ind2~paste(ind1,"_",tis1,".",ind2,"_",tis2,sep="")))
  }
  
  if(mot_cat%in%c("t1g2")){
    fst.obj.tmp<-fst.obj%>%
      mutate(ind1=case_when(ind1=="g1"~"c1",
                            ind1=="g2"~"c2",
                            ind1=="t"~"m1"),
             ind2=case_when(ind2=="g1"~"c1",
                            ind2=="g2"~"c2",
                            ind2=="t"~"m1"))%>%
      mutate(pair=paste(ind1,"_",tis1,".",ind2,"_",tis2,sep=""))
  }
  
  fst.obj.tmp<-fst.obj.tmp%>%
    mutate(branch.type=case_when(pair%in%c("c1_bl.c1_ch","c1_ch.c1_bl","m1_bl.m1_ch","m1_ch.m1_bl","c2_bl.c2_ch","c2_ch.c2_bl")~"bl2ch.div",
                                 pair%in%c("c1_bl.m1_bl","c2_bl.m1_bl")~"m_bl.c_bl",
                                 pair%in%c("c1_bl.m1_ch","c2_bl.m1_ch","c1_ch.m1_bl","c2_ch.m1_bl")~"m_bl.c_ch",
                                 pair%in%c("c1_ch.m1_ch","c2_ch.m1_ch")~"m_ch.c_ch",
                                 pair%in%c("c1_bl.c2_bl")~"c_bl.c_bl",
                                 pair%in%c("c1_ch.c2_ch")~"c_ch.c_ch",
                                 pair%in%c("c1_bl.c2_ch","c1_ch.c2_bl")~"c_bl.c_ch"))
  
  return(fst.obj.tmp)
  
}

#write function to sample 2 kids from each mother (if she has more than 2)
#oldest and youngest of the two to get more divergence and more power
#also reformat dataframe for bls calculations
#this is useful if fst needs to be calculated using the ratio of averages method
oldnyoung<-function(fst.obj){
  mot_cat=unique(fst.obj$mot_cat)
  #if 2 kids, return without changing
  if(mot_cat==c("m1c2")){fst.obj.tmp<-fst.obj}
  #if>2kids, sample 1 and rename them to be c1 and c2 based on age
  if(mot_cat%in%c("m1c3","m1c4","m1c5")){
    children.levels<-setdiff(fst.obj$ind1,"m1")
    # sampled.children<-sort(sample(children.levels,2,replace=F))
    sampled.children<-c(children.levels[1], children.levels[length(children.levels)])
    fst.obj.tmp<-fst.obj[which(fst.obj$ind1%in%c("m1",sampled.children) &
                                 fst.obj$ind2%in%c("m1",sampled.children)),]
    fst.obj.tmp<-fst.obj.tmp%>%
      mutate(ind1=case_when(ind1==sampled.children[1]~"c1",
                            ind1==sampled.children[2]~"c2",
                            !ind1%in%sampled.children~ind1),
             ind2=case_when(ind2==sampled.children[1]~"c1",
                            ind2==sampled.children[2]~"c2",
                            !ind2%in%sampled.children~ind2))%>%
      mutate(pair=paste(ind1,"_",tis1,".",ind2,"_",tis2,sep=""))
  }
  #standardize notation for g1m2 and t1g2
  if(mot_cat%in%c("g1m2")){
    fst.obj.tmp<-fst.obj%>%
      mutate(ind1=case_when(ind1=="m1"~"c1",
                            ind1=="m2"~"c2",
                            ind1=="g1"~"m1"),
             ind2=case_when(ind2=="m1"~"c1",
                            ind2=="m2"~"c2",
                            ind2=="g1"~"m1"))%>%
      mutate(pair=case_when(ind1=="c1"&ind2=="c2"~paste(ind1,"_",tis1,".",ind2,"_",tis2,sep=""),
                            ind1=="m1"&ind2=="c1"~paste(ind2,"_",tis2,".",ind1,"_",tis1,sep=""),
                            ind1=="m1"&ind2=="c2"~paste(ind2,"_",tis2,".",ind1,"_",tis1,sep=""),
                            ind1==ind2~paste(ind1,"_",tis1,".",ind2,"_",tis2,sep="")))
  }
  
  if(mot_cat%in%c("t1g2")){
    fst.obj.tmp<-fst.obj%>%
      mutate(ind1=case_when(ind1=="g1"~"c1",
                            ind1=="g2"~"c2",
                            ind1=="t"~"m1"),
             ind2=case_when(ind2=="g1"~"c1",
                            ind2=="g2"~"c2",
                            ind2=="t"~"m1"))%>%
      mutate(pair=paste(ind1,"_",tis1,".",ind2,"_",tis2,sep=""))
  }
  
  fst.obj.tmp<-fst.obj.tmp%>%
    mutate(branch.type=case_when(pair%in%c("c1_bl.c1_ch","c1_ch.c1_bl","m1_bl.m1_ch","m1_ch.m1_bl","c2_bl.c2_ch","c2_ch.c2_bl")~"bl2ch.div",
                                 pair%in%c("c1_bl.m1_bl","c2_bl.m1_bl")~"m_bl.c_bl",
                                 pair%in%c("c1_bl.m1_ch","c2_bl.m1_ch","c1_ch.m1_bl","c2_ch.m1_bl")~"m_bl.c_ch",
                                 pair%in%c("c1_ch.m1_ch","c2_ch.m1_ch")~"m_ch.c_ch",
                                 pair%in%c("c1_bl.c2_bl")~"c_bl.c_bl",
                                 pair%in%c("c1_ch.c2_ch")~"c_ch.c_ch",
                                 pair%in%c("c1_bl.c2_ch","c1_ch.c2_bl")~"c_bl.c_ch"))
  
  return(fst.obj.tmp)
  
}

"%ni%"<-Negate("%in%")

