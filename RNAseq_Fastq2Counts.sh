
#!/bin/bash
#SBATCH -J RNAseq
#SBATCH -c 8
#SBATCH -p xeon
#SBATCH --qos=general
#SBATCH --mem=100G
#SBATCH -o %x-%j.out


while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo " "
			echo "If the script fails to run to the end, look for the error correct the erreo and rerun the script."
			echo "The script will rerun from the failed step than from the start.   ENJOY.."
			echo " Syntax :"
			echo "		sbatch RNAseq_Fastq2Counts.sh -s <SampleName> -p <Path/to/project directory> -m <human|mouse> "
			echo " "
			echo "OPTIONS: "
			echo "		-s SampleName Please read the the project directory set up below "
			echo "		-p Absolute Path to the project directory "
			echo "		-m human mouse : At present the script supports only these two. For any other species set IndexPath variable in script "
			echo " "
			echo "Directory Set Up:"
			echo "Inside the project directory there is raw_data directory containg sample directories with fastq.gz read data from each sample"
			echo "sample directories should have fastq files"
			echo "ProjectDirectory "
			echo "       |-raw_data "
			echo "               |-samplename1 Directory  "
			echo "               |       |-fastq files "
			echo "               |-samplename2 Directory "
			echo "                       |-fastq files "
			echo ""
			echo "DONOT PUT ALL FASTQ FILES IN SAME DIRECTORY "
			echo " " 
			exit 1 
			;;
                -s)
                        shift
                        if test $# -gt 0; then
                                export sample=$1
                        else
                                echo "no sample specified"
                                exit 1
                        fi
                        shift
                        ;;
                -p)
                        shift
                        if test $# -gt 0; then
                                export result_dir=$1
                        else
                                echo "no Project directory specified"
                                exit 1
                        fi
                        shift
                        ;;
		-m)
                        shift
                        if test $# -gt 0; then
                                export model=$1
                        else
                                echo "human or mouse: not  specified"
                                exit 1
                        fi
                        shift
                        ;; 
                *)
                        break
                        ;;
        esac
done

#echo $sample
#echo $projectDir
#echo $model 
ProjectDir=$result_dir

mkdir -p ${ProjectDir}/tmp
export TMPDIR=${ProjectDir}/tmp

hostname

echo "SampleID :  ${sample}"

#Paths below are specific to Xanadu cluster and can be changed to users need.
# Future development:  Make it more flexible and source through a config file.

if [ $model == human ];then
 	IndexPath="/isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_HISAT2/"
 	IndexPrefix="${IndexPath}/hg19"
elif [ $model == mouse ];then
 	IndexPath="/isg/shared/databases/alignerIndex/animal/mus_musculus/HISAT2/"
 	IndexPrefix="${IndexPath}/Mus_musculus"
else
 	echo " model information not provided.  Set 3rd argument to human or mouse."  
 	echo "IF EITHER OF THEM IS NOT CORRECT then SET THE PATH TO INDEX LOCATION (down in CMD9). "
 	exit 1
fi

if [ ! -d ${ProjectDir}/logs]; then
        mkdir -p ${ProjectDir}/logs
fi

touch ${ProjectDir}/logs/${sample}_OUT
echo "" > ${ProjectDir}/logs/${sample}_OUT

if [ ! -e "${ProjectDir}/logs/${sample}_log_trace" ]; then
        touch ${ProjectDir}/logs/${sample}_log_trace
 	for i in $(seq 1 14); do echo "CMD${i}=FAIL" >> ${ProjectDir}/logs/${sample}_log_trace; done
fi

cd ${ProjectDir}/raw_data/${sample}

OUT=${ProjectDir}/logs/${sample}_OUT
source ${ProjectDir}/logs/${sample}_log_trace

#TEMPLATE
##########################################################################################################
#if [ CMD1 == "OK" ];then
# echo  "Merging of R1 reads was successfully executed in previous run " >> ${OUT}
#else
# $COMMAND
# if [ $? == 0 ]; then
#      sed -i -e 's/CMD1=FAIL/CMD1=OK/g' ${ProjectDir}/logs/${sample}_log_trace
# else
#      echo "CMD1, Mergeing R1, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
# fi
#fi
##########################################################################################################

files_gz=`ls *gz 2>/dev/null | wc -l`

#CMD1
if [ ${CMD1} == "OK" -o  ${files_gz} -eq 0 ];then
  	echo  "UNCOMPRESSSING STEP: Either Uncompressing of fastq files was successfully executed in previous run OR there are no compressed files. " >> ${OUT}
	sed -i -e 's/CMD1=FAIL/CMD1=OK/g' ${ProjectDir}/logs/${sample}_log_trace
else
  	gunzip *.fastq.gz
  	if [ $? == 0 ]; then
       		sed -i -e 's/CMD1=FAIL/CMD1=OK/g' ${ProjectDir}/logs/${sample}_log_trace
  	else
       		echo "CMD1, Uncompressing fastq files, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
   		exit 1
 	fi
fi

NumberOfFiles=`ls *.fastq 2>/dev/null | wc -l`
if [ $NumberOfFiles -eq 8 ];then
  
#CMD2
##########################################################################################################
 	if [ ${CMD2} == "OK" ];then
  		echo  "MERGING R1 FASTQ FILES: Merging of R1 reads was successfully executed in previous run OR NOT REQUIRED " >> ${OUT}
 	else
  #cat ${sample}_S*_L001_R1_001.fastq ${sample}_S*_L002_R1_001.fastq ${sample}_S*_L003_R1_001.fastq ${sample}_S*_L004_R1_001.fastq >> ${sample}_R1.fastq
  		cat `ls *R1*.fastq | sort -` >> ${sample}_R1.fastq
  		if [ $? == 0 ]; then
       			sed -i -e 's/CMD2=FAIL/CMD2=OK/g' ${ProjectDir}/logs/${sample}_log_trace
  		else
       			echo "CMD2, Mergeing R1, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
   			exit 1
  		fi
 	fi
##########################################################################################################

#CMD3
##########################################################################################################
	if [ ${CMD3} == "OK" ];then
 		echo  "MERGING R2 FASTQ FILES: Merging of R2 reads was successfully executed in previous run OR NOT REQUIRED " >> ${OUT}
	else
 #cat ${sample}_S*_L001_R2_001.fastq ${sample}_S*_L002_R2_001.fastq ${sample}_S*_L003_R2_001.fastq ${sample}_S*_L004_R2_001.fastq >> ${sample}_R2.fastq
 	cat `ls *R2*.fastq | sort -` >> ${sample}_R2.fastq
 		if [ $? == 0 ]; then
      			sed -i -e 's/CMD3=FAIL/CMD3=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 		else
      			echo "CMD3, Mergeing R2, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  			exit 1
 		fi
	fi
##########################################################################################################

else
	sed -i -e 's/CMD2=FAIL/CMD2=OK/g' ${ProjectDir}/logs/${sample}_log_trace
	sed -i -e 's/CMD3=FAIL/CMD3=OK/g' ${ProjectDir}/logs/${sample}_log_trace
fi


mkdir -p ${ProjectDir}/merged_fastq

mv ${sample}_R?.fastq ${ProjectDir}/merged_fastq/ 2>/dev/null

cd ${ProjectDir}/merged_fastq

#CMD4
##########################################################################################################
if [ ${CMD4} == "OK" ];then
 	echo  "VALIDATING READS IN PAIRED FILES: File merging Authenticated in previous run " >> ${OUT}
else
 	reads1=`wc -l < ${sample}_R1.fastq`
 	reads2=`wc -l < ${sample}_R2.fastq`
 	total_reads=$(( $reads1 / 4 ))
 	if [ $reads1 == $reads2 ];then
         	echo "File merging Authenticated"
         	echo "Total reads are $total_reads"
 	else
         	echo "There was issue merging the files, The number of reads in both files donot match...!!!!"
  #sed -i -e 's/CMD4=FAIL/CMD4=OK/g' ${ProjectDir}/logs/${sample}_log_trace
         	exit 1
 	fi
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD4=FAIL/CMD4=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD4, Read counts in pair, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################


mkdir -p ${ProjectDir}/fastqc_raw

module load fastqc

#CMD5
##########################################################################################################
if [ ${CMD5} == "OK" ];then
 	echo  "FASTQC, quality check was executed in previous run " >> ${OUT}
else
 	fastqc -t 8 -o ${ProjectDir}/fastqc_raw ${ProjectDir}/merged_fastq/${sample}_R1.fastq ${ProjectDir}/merged_fastq/${sample}_R2.fastq
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD5=FAIL/CMD5=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD5, FastQC, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################

module rm fastqc

mkdir -p ${ProjectDir}/trimmed_reads

module load sickle

#CMD6
##########################################################################################################
if [ ${CMD6} == "OK" ];then
 	echo  "SICKLE : Read trimming was successfully executed in previous run " >> ${OUT}
else
 sickle pe -t sanger \
	-f ${ProjectDir}/merged_fastq/${sample}_R1.fastq \
	-r ${ProjectDir}/merged_fastq/${sample}_R2.fastq \
        -o ${ProjectDir}/trimmed_reads/trimmed_${sample}_R1.fastq \
        -p ${ProjectDir}/trimmed_reads/trimmed_${sample}_R2.fastq \
        -l 45 \
        -q 25 \
        -s ${ProjectDir}/trimmed_reads/singles_${sample}_R1.fastq
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD6=FAIL/CMD6=OK/g' ${ProjectDir}/logs/${sample}_log_trace
	 else
      		echo "CMD6, Read Trimming, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################

cd ${ProjectDir}/trimmed_reads/

#CMD7
##########################################################################################################
if [ ${CMD7} == "OK" ];then
 	echo  "VALIDATING READ PAIRS POST SICKLE: Read counts post Trimming was successfully executed in previous run " >> ${OUT}
else
 	reads_trm1=`wc -l < trimmed_${sample}_R1.fastq`
 	reads_trm2=`wc -l < trimmed_${sample}_R2.fastq`
 	if [ $reads_trm1 == $reads_trm2 ];then
         	echo "Processed sickle files are Authenticated"
 	else
         	echo "Processed sickle output files cannot be authenticated, The number of reads in both files donot match...!!!!"
        exit 1
 	fi
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD7=FAIL/CMD7=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD7, Read counts post Trimming , FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################



mkdir -p ${ProjectDir}/fastqc_processed

module rm sickle

module load fastqc

#CMD8
##########################################################################################################
if [ ${CMD8} == "OK" ];then
 	echo  "FastQC post trimming was successfully executed in previous run " >> ${OUT}
else
 	fastqc -t 8 -o ${ProjectDir}/fastqc_processed ${ProjectDir}/trimmed_reads/trimmed_${sample}_R1.fastq ${ProjectDir}/trimmed_reads/trimmed_${sample}_R2.fastq
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD8=FAIL/CMD8=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD8, FastQC post trimming , FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################

module rm fastqc

mkdir -p ${ProjectDir}/mapping
mkdir -p ${ProjectDir}/mapping/tmp_${sample}

module load hisat2/2.1.0

#CMD9
##########################################################################################################
if [ ${CMD9} == "OK" ];then
 	echo  "HiSAT2 mapping of reads was successfully executed in previous run " >> ${OUT}
else
	hisat2 -p 8 --known-splicesite-infile ${IndexPath}/known_splice_sites.txt \
	-x ${IndexPrefix} \
	-1 ${ProjectDir}/trimmed_reads/trimmed_${sample}_R1.fastq \
	-2 ${ProjectDir}/trimmed_reads/trimmed_${sample}_R2.fastq \
	-S ${ProjectDir}/mapping/tmp_${sample}/${sample}.sam

 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD9=FAIL/CMD9=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
		echo "CMD9, HISAT2 mapping, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
		echo "Note: If (ERR): hisat2-align died with signal 11 (SEGV) Is the fail error." >>${ProjectDir}/logs/${sample}_OUT
		echo "Check hostnode, mentioned at the top of SLURM OUTPUT, it will be either shangrila or xandau" >>${ProjectDir}/logs/${sample}_OUT
	 	echo "If it is shangrila, SIMPLY RERUN THE CODE, There is issue with memory clearance in shangrila nodes" >>${ProjectDir}/logs/${sample}_OUT
  	 	exit 1
 	fi
fi
##########################################################################################################

cd ${ProjectDir}/mapping/tmp_${sample}

module load samtools/1.7

#CMD10
##########################################################################################################
if [ ${CMD10} == "OK" ];then
 	echo  "SAM -> BAM conversion  was successfully executed in previous run " >> ${OUT}
else
 	samtools view -@ 8 -bhS ${sample}.sam -o ${sample}_mapped.bam
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD10=FAIL/CMD10=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD10, Samtool SAM -> BAM conversion, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################

#CMD11
##########################################################################################################
if [ ${CMD11} == "OK" ];then
 	echo  "BAM file sorting (Before PCR duplicate removal) was successfully executed in previous run " >> ${OUT}
else
 	samtools sort -@ 8 ${sample}_mapped.bam -o ${sample}_mapped_sort.bam
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD11=FAIL/CMD11=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD11, Samtool sorting of .bam files, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################


module load java
#CMD12
##########################################################################################################
if [ ${CMD12} == "OK" ];then
 	echo  "PCR Duplicate removal was successfully executed in previous run " >> ${OUT}
else
 	java -jar /isg/shared/apps/picard/picard-tools-2.2.1/picard.jar MarkDuplicates INPUT= ${ProjectDir}/mapping/tmp_${sample}/${sample}_mapped_sort.bam OUTPUT= ${ProjectDir}/mapping/tmp_${sample}/${sample}_dup_removed.bam METRICS_FILE= ${sample}_dup_removed_metrics REMOVE_DUPLICATES= true CREATE_INDEX= true
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD12=FAIL/CMD12=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD12, Mergeing R1, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################
#CMD13
##########################################################################################################
if [ ${CMD13} == "OK" ];then
 	echo  "BAM file sorting, post PCR duplicate removal, was successfully executed in previous run " >> ${OUT}
else
 	samtools sort -@ 8 ${sample}_dup_removed.bam -o ${sample}_dup_removed_sort.bam
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD13=FAIL/CMD13=OK/g' ${ProjectDir}/logs/${sample}_log_trace
 	else
      		echo "CMD13, BAM file sorting post PCR duplicate removal, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
  		exit 1
 	fi
fi
##########################################################################################################


mv ${sample}_dup_removed_sort.bam ${ProjectDir}/mapping/

module rm samtools/1.7
module rm hisat2/2.1.0

mkdir -p ${ProjectDir}/counts

cd ${ProjectDir}/counts

module load htseq/0.11.0
#CMD14
##########################################################################################################
if [ ${CMD14} == "OK" ];then
 	echo  "Merging of R1 reads was successfully executed in previous run " >> ${OUT}
else

 	htseq-count -s no -r pos -t exon -i gene_id -f bam ${ProjectDir}/mapping/${sample}_dup_removed_sort.bam /isg/shared/databases/alignerIndex/animal/mus_musculus/Mus_musculus.GRCm38.93.gtf > ${ProjectDir}/counts/${sample}.counts
 	if [ $? == 0 ]; then
      		sed -i -e 's/CMD14=FAIL/CMD14=OK/g' ${ProjectDir}/logs/${sample}_log_trace
	else
      		echo "CMD14, Mergeing R1, FAILED CHECK THE SLURM OUTPUT FOR DETAILS" >>${ProjectDir}/logs/${sample}_OUT
	  	exit 1
 	fi
fi
##########################################################################################################

module rm htseq/0.11.0

echo " =====================================================" >> ${ProjectDir}/logs/${sample}_OUT
echo " " >> ${ProjectDir}/logs/${sample}_OUT
echo " ANALYSIS COMPLETED" >> ${ProjectDir}/logs/${sample}_OUT

##########################################################################################################

echo "" >> ${ProjectDir}/logs/${sample}_OUT
echo "" >> ${ProjectDir}/logs/${sample}_OUT


