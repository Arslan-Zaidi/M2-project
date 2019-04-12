#!/usr/bin/env python

import pandas as pd
import os
import argparse
from scipy.stats import poisson

parser = argparse.ArgumentParser()
parser.add_argument('-c',dest='count_file',required=True,help='path to count file')
parser.add_argument('-f',dest='maf',default=0.01,help='minor allele frequency threshold. default=0.01',type=float)
parser.add_argument('-d',dest='depth',default=1000,help='minimum depth of dequencing. default = 1000',type=int)
parser.add_argument('-o',dest='output',required=True,help='output filename')

args = parser.parse_args()

print "reading in count file"


df=pd.read_table(args.count_file,names=["fqid","reference","position","A","C","G","T","a","c","g","t","cvrg","nalleles","major","minor","maf"],dtype={"fqid":str, "reference":str,"position":int,"A":int,"C":int,"G":int,"T":int,"a":int,"c":int,"g":int,"t":int,"cvrg":int,"nalleles":int,"major":str,"minor":str,"maf":float},header=None)
df['tr']=args.count_file.split(".")[0]


##cleaning up count file.
#df.columns=["fqid","reference","position","A","C","G","T","a","c","g","t","cvrg","nalleles","major","minor","maf","tr"]

#filtering out clones and un-indexed reads
df2=df[-df['fqid'].str.startswith('Z1')]
df2=df2[-df2['fqid'].str.startswith('Unde')]
df2=df2[-df2['fqid'].str.startswith('z1')]



##establishing filters for heterolplasmy detection
#setting up to filter hq_sites

#1. masking 'bad' regions

#low complexity regions: chrM:66-71, chrM:303-311, chrM:514-523, chrM:12418-12425, chrM:16183-16193, chrM:3105-3109
#50bp downstream and 25bp upstream of primer binding sites: 
#LRA_F_primer: 2817..2868, LRA_R_primer: 11520..11570
#LRB_F_primer: 3320..3370, LRB_R_primer: 10796..10846
mask = [(66,71),(303,311),(514,523),(12418,12425),(16183,16193),
      (3105,3109),(2792,2868),(3320,3370),(10796,10846),(11520,11570)]


maskRegions = list()
for start,end in mask:
    maskRegions+=range(start,end+1)

    
#2. define function to calculate maf bias and to filter on maf bias
def strand_stats(x, mafThreshold=args.maf):
    falleles = ['A','C','G','T']
    ralleles = ['a','c','g','t']
    sample,position,major,minor,coverage,maf = x[['fqid','position','major','minor','coverage','maf']]
    fcounts = x[falleles]
    rcounts = x[ralleles]
    if minor!='.':
        index_major = falleles.index(major)
        index_minor = falleles.index(minor)

        fcount_minor = float(fcounts[index_minor])
        ftotal = fcount_minor + fcounts[index_major]
        
        rcount_minor = float(rcounts[index_minor])
        rtotal = rcount_minor + rcounts[index_major]
        
        minor_total = float(fcount_minor + rcount_minor)
        site_total = ftotal + rtotal

        try:
            strandBias = abs( (fcount_minor/ftotal) - (rcount_minor/rtotal) ) / (minor_total/site_total)
        except:
            strandBias = np.nan
            
        try:
            maf_frwd = fcount_minor/sum(fcounts)
        except:
            maf_frwd = np.nan
        try:
            maf_rvrs = rcount_minor/sum(rcounts)
        except:
            maf_rvrs = np.nan
            
        if (maf_frwd>=mafThreshold) and (maf_rvrs>=mafThreshold):
            mafBalance = 1
        else:
            mafBalance = 0
    else:
        strandBias = float(2)
        mafBalance = 0

    return pd.Series([strandBias,mafBalance])



#3. define function to calculate probability of each heteroplasmy based on a poisson distribution

def poisson_pval(current_df,sample):
    alleles = ['A','C','G','T','a','c','g','t']
    
    sample_counts = list(current_df.loc[current_df['fqid']==sample, alleles].iloc[0,:])
    others_counts = list(current_df.loc[current_df['fqid']!=sample, alleles].apply(sum,axis=0))
    sample_coverage = sum(sample_counts)
    
    observed_error = (sum(others_counts) - max(others_counts))/float(sum(others_counts))
    sample_nonMajor_counts = int(sample_coverage - max(sample_counts))
    
    pvalue = poisson.pmf(sample_nonMajor_counts, observed_error*sample_coverage)
    
    return pvalue

 ##applying filters to df to generate hq_sites

print "Filtering count table on minor allele frequency and depth"
#filter on maf (>=1%) and coverage (>=2000x)
hq= df2[(df2.maf>=args.maf) & (df2.cvrg>=args.depth) & ~df2.position.isin(maskRegions)]


print "Filtering on strand bias and maf balance"

#filter on strand bias and maf balance
biasCols = hq.apply(strand_stats,axis=1)
biasCols.columns = ["strandBias","mafBalance"]
hq = pd.concat([hq,biasCols],axis=1)
hq= hq[(hq.strandBias<=1) & (hq.mafBalance==1) ]

print "Excluding variants on Poisson p-values"
# filter on possion probability

poisson_pvalues = []

for sample,position in hq[["fqid","position"]].itertuples(index=False):
    poisson_pvalues.append(poisson_pval(df2[df2['position']==position],sample))
    
hq["poisson"] = poisson_pvalues
hq = hq[hq.poisson<=0.05]

print "Writing table of high quality heteroplasmies to :",args.output 

hq.to_csv(args.output,index=False,header=True,sep="\t")