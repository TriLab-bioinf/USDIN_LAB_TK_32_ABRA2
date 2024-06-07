### USDIN_LAB_TK_32_ABRA2
Updated SVN annotation pipeline based on BWA, ABRA2, GATK and snpEff

**SVN_calling_pipeline_bwa_abra2_gatk**
This folder contains the scripts required for the identification of SNPs and INDELs in one or more samples. To run the pipeline you need first to configure the *config.cfg* file with the following information:

```
# path to reads directory
READ=/gpfs/gsfs12/users/lorenziha/KAREN_USDIN/TK_119/USDIN_LAB_TK_32/comparative_analysis/READS/

# path to adapters fasta file
ADAPTERS=/gpfs/gsfs12/users/lorenziha/KAREN_USDIN/TK_119/USDIN_LAB_TK_32/adapters.fa

# Path to reference genome
REFERENCE=/fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa

# Merge reads before mapping? [true (default) or false]
MERGE=false
```
Then you need to edit the line below from the *snpEff.config* file to indicate the path where to create the snpEff annotation file.
**IMPORTANT:** The *snpEff.config* file has to be in the same directory as the *bruce_pipeline_2.sh* script.
```
data.dir = /gpfs/gsfs12/users/lorenziha/KAREN_USDIN/TK_119/USDIN_LAB_TK_32/snpeff_db
```
Next edit the *prefix.txt* file by adding the prefix of the samples to be processed, one prefix per line.
**IMPORTANT:** It is expected that fastq file names follow the following naming convention:
```
<prefix>.R1.fastq.gz
<prefix>.R2.fastq.gz
```
Finally, to run the pipeline in Biowulf type the following command:
```
sbatch --time=24:00:00 bruce_pipeline_2.sh prefix.txt config.cfg
```

**Rstudio.server.scripts** folder contains the scripts to run ExonDepth analysis (identification of CNVs) using R in *Rstudio.server* from Biowulf.

As implemented right now, the analysis has to be run interactively from Rstudio.server in Biowulf using the following R script: *run_exomeDepth_server.R*. The R script *aux1-exomedepth.R* contains the actual ExomeDepth function and generated a csv file summarizing the results and two ideogram plots in a folder named *07-exomeDepth*.

Before running you will need to edit the file *control_bams.txt* with the path to each of the control bam files. Also, you will need to adjust the values of the *case_files* and *input_bam_file* variables to point to the correct sample bam file. It is expected that sample bam files are located in folder *04-abra2* and follow the following name convention:
```
<SAMPLE PREFIX>_realigned.bam
```






