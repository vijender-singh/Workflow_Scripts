#!/bin/bash
#SBATCH -J RNAseq
#SBATCH -c 8
#SBATCH -p general
#SBATCH --qos=general
#SBATCH --mem=100G
#SBATCH -o RNAseq-%x-%j.out

##MAKESURE raw_data DIRECTORY IS INSIDE $ProjectDir AND IT HAS SAMPLE DIRECTORIES CONTAINNING FASTQ FILES INSIDE THEM.
## ALSO WITHIN ProjectDir MAKE A script DIRECTORY WITH known_splice_sites.txt EXTRACTED FROM GFF OR GTF FILE

ProjectDir=$1

mkdir -p ${ProjectDir}/tmp

export TMPDIR=${ProjectDir}/tmp

sample="$1"


cd ${ProjectDir}/raw_data/${sample}

gunzip *.fastq.gz

hostname

echo ${sample}

cat ${sample}_S*_L001_R1_001.fastq ${sample}_S*_L002_R1_001.fastq ${sample}_S*_L003_R1_001.fastq ${sample}_S*_L004_R1_001.fastq >> ${sample}_R1.fastq
cat ${sample}_S*_L001_R2_001.fastq ${sample}_S*_L002_R2_001.fastq ${sample}_S*_L003_R2_001.fastq ${sample}_S*_L004_R2_001.fastq >> ${sample}_R2.fastq

mkdir -p ${ProjectDir}/merged_fastq

mv ${sample}_R?.fastq ${ProjectDir}/merged_fastq/

cd ${ProjectDir}/merged_fastq

reads1=`wc -l < ${sample}_R1.fastq`
reads2=`wc -l < ${sample}_R2.fastq`
total_reads=$(( $reads1 / 4 ))
if [ $reads1 == $reads2 ];then
        echo "File merging Authenticated"
        echo "Total reads are $total_reads"
else
        echo "There was issue merging the files, The number of reads in both files donot match...!!!!"
        exit 1
fi

mkdir -p ${ProjectDir}/fastqc_raw

module load fastqc

fastqc -t 8 -o ${ProjectDir}/fastqc_raw ${ProjectDir}/merged_fastq/${sample}_R1.fastq ${ProjectDir}/merged_fastq/${sample}_R2.fastq

module rm fastqc

mkdir -p ${ProjectDir}/trimmed_reads

module load sickle
sickle pe -t sanger \
        -f ${ProjectDir}/merged_fastq/${sample}_R1.fastq \       
	-r ${ProjectDir}/merged_fastq/${sample}_R2.fastq \
        -o ${ProjectDir}/trimmed_reads/trimmed_${sample}_R1.fastq \
        -p ${ProjectDir}/trimmed_reads/trimmed_${sample}_R2.fastq \
        -l 45 \
        -q 25 \
        -s ${ProjectDir}/trimmed_reads/singles_${sample}_R1.fastq

cd ${ProjectDir}/trimmed_reads/

reads_trm1=`wc -l < trimmed_${sample}_R1.fastq`
reads_trm2=`wc -l < trimmed_${sample}_R2.fastq`

if [ $reads_trm1 == $reads_trm2 ];then
        echo "Processed sickle files are Authenticated"
else
        echo "Processed sickle output files cannot be authenticated, The number of reads in both files donot match...!!!!"
        exit 1
fi


mkdir -p ${ProjectDir}/fastqc_processed

module rm sickle

module load fastqc

fastqc -t 8 -o ${ProjectDir}/fastqc_processed ${ProjectDir}/trimmed_reads/trimmed_${sample}_R1.fastq ${ProjectDir}/trimmed_reads/trimmed_${sample}_R2.fastq

module rm fastqc

mkdir -p ${ProjectDir}/mapping
mkdir -p ${ProjectDir}/mapping/tmp_${sample}

module load hisat2/2.1.0

hisat2 -p 8 --known-splicesite-infile ${ProjectDir}/scripts/known_splice_sites.txt \
        -x /isg/shared/databases/alignerIndex/animal/mus_musculus/HISAT2/Mus_musculus \ # CORRECT THIS PATH IF NOT RIGHT
        -1 ${ProjectDir}/trimmed_reads/trimmed_${sample}_R1.fastq \
        -2 ${ProjectDir}/trimmed_reads/trimmed_${sample}_R2.fastq \
        -S ${ProjectDir}/mapping/tmp_${sample}/${sample}.sam

cd ${ProjectDir}/mapping/tmp_${sample}

module load samtools/1.7
samtools view -@ 8 -bhS ${sample}.sam -o ${sample}_mapped.bam
samtools sort -@ 8 ${sample}_mapped.bam -o ${sample}_mapped_sort.bam

module load java

java -jar /isg/shared/apps/picard/picard-tools-2.2.1/picard.jar MarkDuplicates INPUT= ${ProjectDir}/mapping/tmp_${sample}/${sample}_mapped_sort.bam OUTPUT= ${ProjectDir}/mapping/tmp_${sample}/${sample}_dup_removed.bam METRICS_FILE= ${sample}_dup_removed_metrics REMOVE_DUPLICATES= true CREATE_INDEX= true

samtools sort -@ 8 ${sample}_dup_removed.bam -o ${sample}_dup_removed_sort.bam

mv ${sample}_dup_removed_sort.bam ${ProjectDir}/mapping/

module rm samtools/1.7
module rm hisat2/2.1.0

mkdir -p ${ProjectDir}/counts

cd ${ProjectDir}/counts

module load htseq/0.11.0

htseq-count -s no -r pos -t exon -i gene_id -f bam ${ProjectDir}/mapping/${sample}_dup_removed_sort.bam /isg/shared/databases/alignerIndex/animal/mus_musculus/Mus_musculus.GRCm38.93.gtf > ${ProjectDir}/counts/${sample}.counts

module rm htseq/0.11.0
