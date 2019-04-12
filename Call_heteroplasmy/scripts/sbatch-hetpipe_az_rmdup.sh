#!/bin/bash  

#SBATCH -C new
##SBATCH --cpus-per-task=8
#SBATCH -J TR2.rd
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=saz5078@psu.edu
#SBATCH -e slurm-%j.err
#SBATCH -D /nfs/brubeck.bx.psu.edu/scratch6/arslan/slrm/M2

source  /nfs/brubeck.bx.psu.edu/scratch4/arslan/.bash_profile
source activate cenv

cd $SLURM_SUBMIT_DIR

wd=${1}
cd $wd
run=${2}
subject=${3}


###############################################################################

reference="/nfs/brubeck.bx.psu.edu/scratch4/arslan/mtproj/mtproject/refgen/mtprojectref"

f1=${run}/${subject}_L001_R1_001.fastq.gz
f2=${run}/${subject}_L001_R2_001.fastq.gz

mkdir -p ${subject}/fq
mkdir -p ${subject}/fil.mapped

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${subject}/fq/${subject}_1.fq -p ${subject}/fq/${subject}_2.fq ${f1} ${f2}

bwa mem \
    -t 8 \
    -M   \
    -R "@RG\tID:${subject}\tLB:amelia\tPL:ILLUMINA\tSM:${subject}"\
    ${reference}      \
    ${subject}/fq/${subject}_1.fq             \
    ${subject}/fq/${subject}_2.fq 2>/dev/null|\
    samtools view -Sb - > ${subject}/fil.mapped/${subject}.mapped.bam 2>/dev/null

###############################################################################

cd ${subject}/fil.mapped
file=${subject}.mapped.bam
name=${subject}


#regex="[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9-]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*."
regex="null"
ref="/nfs/brubeck.bx.psu.edu/scratch4/arslan/rcrs/rcrs.fasta"
#picardTools="/nfs/brubeck.bx.psu.edu/scratch4/arslan/bin"

###############################################################################

echo "----determine read length started..."

readLength=`samtools view ${file} |head -100|awk '{sum+=length($10)} END{print sum/NR}'`

echo "----determine read length finished..."
###############################################################################
#picard=TMP_DIR=tmp VALIDATION_STRINGENCY=SILENT VERBOSITY=ERROR QUIET=true
echo "----sort by coordinate started..."
picard SortSam \
    ${picard}                        \
    I=${file}                        \
    O=srt.${file}                    \
    SO=coordinate
echo "----sort by coordinate finished..."
###############################################################################

echo "----markDup started..."
picard MarkDuplicates \
    ${picard}                               \
    I=srt.${file}                           \
    REMOVE_DUPLICATES=true                  \
    O=marked.srt.rd.${file}                    \
    METRICS_FILE=srt.${file}.metrics        \
    READ_NAME_REGEX=${regex}
echo "----markDup finished..."
###############################################################################

echo "----Selection of proper started..."
bamtools filter                   \
    -region chrM                  \
    -isPaired true                \
    -isProperPair true            \
    -isMateMapped true            \
    -isPrimaryAlignment true      \
    -in marked.srt.rd.${file}        \
    -out proper.marked.srt.rd.${file}
echo "----Selection of proper finished..."
###############################################################################

echo "----sort by name started..."
picard SortSam  \
    ${picard}                         \
    I=proper.marked.srt.rd.${file}       \
    O=query.proper.marked.srt.rd.${file} \
    SO=queryname
echo "----sort by name finished..."
###############################################################################

echo "----dechim and rlen started..."
rm_chim_in_pair.py query.proper.marked.srt.rd.${file} ${readLength}
echo "----dechim and rlen finished..."
###############################################################################

echo "----realignment started..."
samtools view -b dechim.rlen.query.proper.marked.srt.rd.${file}|               \
bamleftalign -f ${ref}|                                                     \
samtools view -b - > realigned.dechim.rlen.query.proper.marked.srt.rd.${file}
echo "----realignment finished..."
###############################################################################

echo "----sort by coordinate started..."
picard SortSam                   \
    ${picard}                                                   \
    I=realigned.dechim.rlen.query.proper.marked.srt.rd.${file}     \
    O=srt.realigned.dechim.rlen.query.proper.marked.srt.rd.${file} \
    SO=coordinate
echo "----sort by coordinate finished..."
###############################################################################

echo "----index bam started..."
samtools index srt.realigned.dechim.rlen.query.proper.marked.srt.rd.${file}
echo "----index bam finished..."
###############################################################################

echo "----major sequence started..."
#get_major_from_bam.py srt.realigned.dechim.rlen.query.proper.marked.srt.${file}
name2="srt.realigned.dechim.rlen.query.proper.marked.srt.rd.${file}"
angsd              \
    -i ${name2}    \
    -doHetPlas 2   \
    -out ${name2}  \
    -nThreads 4    \
    -minQ 30       \
    -minMapQ 20    \
    -r chrM        \
    -nLines 10000  \
    -doCounts 1    \
    -dumpCounts 3  \
    -doQsDist 1    \
    -howOften 10

mitoMajorFromThorGL_az.py ${name2}.hetGL
echo "----major sequence finished..."
###############################################################################

echo "----recalculate NM started..."
samtools fillmd               \
    -b ${name2}               \
    ${name2}.hetGL.major.fa   \
    2>/dev/null > md.${name2}
echo "----recalculate NM finished..."
###############################################################################

echo "----sort by name 2 started..."
picard SortSam \
    ${picard}                        \
    I=md.${name2}                    \
    O=query.md.${name2}              \
    SO=queryname
echo "----sort by name 2 finished..."
###############################################################################

echo "----NM selection started..."
nm-ratio.select.py query.md.${name2}
echo "----NM selection finished..."
###############################################################################

mv nm-ratio.query.md.${name2} ${name}.rd.nvcReady.bam
###############################################################################

echo "----sort by coordinate started..."
picard SortSam \
    ${picard}                        \
    I=${name}.rd.nvcReady.bam           \
    O=srt.${name}.rd.nvcReady.bam       \
    SO=coordinate
echo "----sort by coordinate finished..."
###############################################################################

echo "----cleaning started..."
mkdir -p nvcReady
mv srt.${name}.rd.nvcReady.bam nvcReady
mkdir -p fasta
mv *.fa fasta/
rm -fr *srt.${file}.* *.nvcReady.*
echo "----cleaning finished..."
###############################################################################

echo "----NVC started..."
cd nvcReady
samtools index srt.${name}.rd.nvcReady.bam

naive_variant_caller.py             \
    -b srt.${name}.rd.nvcReady.bam     \
    -i srt.${name}.rd.nvcReady.bam.bai \
    -o ${name}.rd.vcf                  \
    -r ${ref}                       \
    -s -q 30 -m 20                  \
    -t uint32                       \
    --region chrM &>/dev/null
echo "----NVC finished..."
###############################################################################

echo "----Allele counts started..."
allele-counts.py        \
    -i ${name}.rd.vcf      \
    -o ${name}.rd.counts   \
    -n -s
echo "----Allele counts finished..."
