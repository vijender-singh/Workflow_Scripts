## Helpfile for RNAseq_Fatq2counts.sh

### Script Usage (Xanadu Cluster)

help file : `sh RNAseq_Fastq2Counts.sh -h`

**Syntax**<br/>
```sh
          sbatch RNAseq_Fastq2Counts.sh -s <SampleName> -p <Path/to/project directory> -m <human|mouse>
```
<br/>

*If the script fails due to errors, please check the cluster error/output files and also ProjectDirectory/logs/sample_log_trace file for errors*.<br/>
*After correcting for error resubmit the script and it should run from the step where it failed then from the start .* <br/>
<br/>

**OPTIONS**<br/>
                -s SampleName Please read the the project directory set up below
                -p Absolute Path to the project directory
                -m human or mouse : At presesent the script supports only these two. For any other species set IndexPath variable in script


### Script Purpose
The script takes fastq files and outputs the counts file.  The steps are listed below<br/>
        1. Merging fastq files. Only executes if the **R1** read files are split into 4 based on lanes (L001, L002, L003 and L004)<br/>
        2. `FASTQC` :  Runs fastqc for quality check<br/>
        3. `Sickle` : Trimms the reads based on quality and length of reads,[Deafult q = 25  and length (-l)=45.<br/>
        4. `FASTQC` : Quality check, Post sickle Processing<br/>
        5. `HISAT2` : Mapping reads to reference.<br/>
        6. `htseq-count` : Count reads/fragments corresponding to features in GTF file.<br/>

It also perform some additional steps like,<br/>
        1. Ensures that files generated following a process are correct.<br/>
        2. SAM -> BAM conversions and sorting of BAM files<br/>


### Directory structure required by the script

ProjectDirectory is the main directory, inside which the script expects a `raw_data` directory.  The `raw_data` should have `sample` directories containing fastq files.<br/>
<br/>
Note: fastq files for each samples should be ina different directory. Names of `sample` directories will be included in output files.  Please name them appropriately so that the results are  distingushable from each other.<br/>


### Output

It will create a number of directories during the process containing the output from each step.<br/>

The `.counts` file generated in the `htseq-count` process will be present in `counts ` directory.<br/>

### Resource Header

The first part of the script is the resource request header.  The resources listed below are specific to **Xanadu** cluster hosted at **University of Connecticut**.  Depending on the HPC resource available to you and the scheduler it runs on, you may have to change these parameters.

```sh
#!/bin/bash
#SBATCH -J RNAseq
#SBATCH -c 8
#SBATCH -p xeon
#SBATCH --qos=general
#SBATCH --mem=100G
#SBATCH -o %x-%j.out
```

if you have any questions or suggestions please email : **vijender.singh@uconn.edu**
