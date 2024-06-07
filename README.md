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





