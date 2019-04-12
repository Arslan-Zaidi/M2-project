# Mitochondrial heteroplasmy analysis

The data and code in this repository correspond to the analysis carried out in the following paper. The scripts are heavily commented and Rmarkdown files are included where necessary to walk the reader through most of the analyses that are part of the paper. If you find any bugs, errors, broken links, feel free to contact me or submit a pull request. Please cite the following paper if you decide to use these resources. 

Paper citation: Zaidi, Wilton ...

## Description of folders:

- **Data_files**: Folder contains all the data files needed to run all the analyses. A readme (.md) file is included in the folder that describes what each data file contains and what the different columns mean.

- **Call_heteroplasmy**: Folder contains scripts, which were used to process FASTQ files (paired-end reads), align them, and generate high-confidence heteroplasmy calls. These files were being used with files on a remote server so they cannot be used at present directly on user files, though you could edit them to suit your needs pretty easily. I am currently working to make them more user-friendly.

- **Filtering_heteroplasmies**: Filtering of heteroplasmies based on various criteria as listed in Fig. S3 from the manuscript.

- **Mutation_spectrum**: Code for the anaysis of the mutation spectrum and comparison of heteroplasmy frequency for deleterious vs neutral mutations in human tissues. RMD and rendered .html files are included which show how the analyses were done.

- **Denovo_mutation_analysis**: Heuristic identification of germline and somatic denovo mutations, and analysis of the effect of age on accumulation of these mutations. RMD and rendered .html files are included which show how the analyses were done.

- **Validations**: Analysis of experimental validation of high-frequency heteroplasmies (MAF>10%) using Sanger sequencing and putative de novo germline mutations with droplet digital PCR.

- **BLS**: Code for all the analyses associated with Branch Length Statistics. The source code file "abs_functions_03152019.R" also contains functions to estimate pairwise Fst from allele frequency/count data. Different estimators are provided (e.g. Weir and Cockerham, Hudson, Reynold's). RMD and rendered .html files are included which show how the analyses were done.
