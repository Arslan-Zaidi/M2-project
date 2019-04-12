# Description of data files

#### Author: AAZaidi

## famfile_cleared_conservative_09272018.txt

This file contains information on all the samples (with and without heteroplasmies) analyzed. The column names are pretty straightforward for the most part. In cases where I thought more information was needed, I've provided some details below:

-run: Sequencing run (TRx refers to sequence data generated for this paper, while XT refers to data generated for Rebolledo-Jaramillo et al. 2014).

-fqid: FASTQ ID. This is also the prefix used to store the sequences on the short read archive (SRA).

-removed: refers to sequences that were filtered out. Blank column

-FID: Family ID

-Mother id: Individual's mother's ID

-Level: Level in the pedigree (t=great grandmother, g1=grandmother, etc.)

-spike_in: Which spike-in was added. phiX, puc18, or none

-fam_str: short for family structure. refers to the type of pedigree the individual comes from (e.g. 0-0-1-2 refers to a family with one mother and two children whereas 1-1-2-2 refers to a family with a great-grandmother, her daughter, two grand-daughters, and two great-grand-daughters)

-fam_cat: another way to represent family structure (tg1m2c2=1-1-2-2)

-mot_cat: Each mother was also assigned to a specific category, which was helpful for the BLS analysis. For example, m1c2 refers to a mother with two children

-age_collection: age of the individual at collection

-age_birth: age of the mother when she gave birth to this individual

-twin_info: If an individual is a twin, the individual_id of the twin will be listed here.

======

## hq_cleaned_conservative_09272018.txt

This file contains the heteroplasmy data (nrows=668 for 668 heteroplasmies + header).

======

## hqcounts_cleaned_conservative_09272018.txt
## hq_counts_adj_frequency_09272018.txt

These files contain the nucleotide counts and frequencies for heteroplasmies from all family members of the 'proband'. The difference b/w the first and second files is that the second file contains an additional column: adj.f. This column contains frequencies for the same allele from all family members whereas the maf column only contains the frequency of the minor allele from that individual, which may not be consistent across family members. 

======

## hq_adjf_ancestral_syn_pathogenic_03222019.txt

This file contains additional annotations for heteroplasmies (e.g. gene, synonymous/non-synonymous, pathogenicity etc.)









