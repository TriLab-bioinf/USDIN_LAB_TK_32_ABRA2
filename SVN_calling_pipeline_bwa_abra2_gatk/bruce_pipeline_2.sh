#!/usr/bin/bash
#SBATCH --mem=64g --cpus-per-task=32 --gres=lscratch:10 --time=1-00:00:00 -J WES_pipeline
#From http://seqanswers.com/forums/showthread.php?t=42776

set -e  

# Command example:
# sbatch --time=24:00:00 bruce_pipeline_2.sh <prefix file [one prefix per line]> config.cfg

PREFIX_FILE=$1
CONFIG=$2


if [[ ! -e ${CONFIG} ]]; then
    CONFIG=config.cfg
fi        

source ${CONFIG}

# Check if TMPDIR exists

if [[ ! -d ${TMPDIR} ]]; then
    mkdir ${TMPDIR}
fi        

module load bbtools/38.96 samtools/1.19 GATK/4.5.0.0 picard/3.1.0 bwa bcftools/1.19 snpEff/5.2

#####################################
# Trim adapters from paired end reads
#####################################

if [[ ! -d "01-trimming" ]]; then
    mkdir "01-trimming"
fi


while read SAMPLE
do
    if [[ ! -e FLAG_${SAMPLE}_TRIM ]]; then
    time bbtools bbduk in1=${READ}/${SAMPLE}.R1.fastq.gz in2=${READ}/${SAMPLE}.R2.fastq.gz \
        ref=${ADAPTERS} \
        out1=./01-trimming/${SAMPLE}_R1_trimmed.fastq.gz \
        out2=./01-trimming/${SAMPLE}_R2_trimmed.fastq.gz \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        -Xmx1g > ./01-trimming/${SAMPLE}_adapter_trim.log 2>&1  
    touch FLAG_${SAMPLE}_TRIM        
    fi
 
    # Quality trimming

    # bbduk.sh -Xmx1g in=read1.fq in=read2.fq out1=clean1.fq out2=clean2.fq qtrim=rl trimq=10

    if [[ ! -e FLAG_${SAMPLE}_QC ]]; then
        time bbtools bbduk in1=./01-trimming/${SAMPLE}_R1_trimmed.fastq.gz in2=./01-trimming/${SAMPLE}_R2_trimmed.fastq.gz \
                out1=./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz \
                out2=./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq.gz \
                qtrim=rl trimq=20 -Xmx1g > ./01-trimming/${SAMPLE}_quality_trim.log 2>&1
        touch FLAG_${SAMPLE}_QC
    fi

    #############################
    # Merge PE reads with bbmerge
    #############################

    if [[ ! -d "02-merge" ]]; then
        mkdir "02-merge"
    fi

    # bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

    if [[ ! -e  FLAG_${SAMPLE}_MERGE ]]; then 
        time bbtools bbmerge in1=./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz in2=./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq.gz \
                out=./02-merge/${SAMPLE}_merged.fastq.gz \
                outu1=./02-merge/${SAMPLE}_R1_unmerged.fastq.gz \
                outu2=./02-merge/${SAMPLE}_R2_unmerged.fastq.gz \
                strict=t > ./02-merge/${SAMPLE}_merge.log 2>&1
        touch FLAG_${SAMPLE}_MERGE 
    fi

    REF=`basename ${REFERENCE}`
    REF_BASENAME=`echo ${REF%.f*a}`

    ############################################################################################
    # Ensure that reference fasta file is in current directory, otherwise create a symbolic link
    ############################################################################################

    if [[ ! -e ${REF} ]]; then
        ln -s ${REFERENCE} ${REF}
    fi          

    #########################################
    # Ensure that reference dict is present 
    #########################################

    if [[ ! -e  ${REF_BASENAME}.dict ]]; then
        samtools dict ${REFERENCE} -o ${REF_BASENAME}.dict       
    fi

    ##############################################
    # Ensure that reference index files is present
    ##############################################

    if [[ ! -e  ${REF}.fai ]]; then
        samtools faidx ${REF} -o ${REF}.fai
    fi        

    ##################################
    # Make mm10 reference genome file
    ##################################

    if [[ ! -e ${REF}.sa  ]]; then
        echo - Making BWA index files
        bwa index ${REF}
    fi            


    #############################################################################
    # Mapping (no need to specify genome as mm10 the only one available)
    #############################################################################

    if [[ ! -d "03-mapping" ]]; then
        mkdir "03-mapping"
    fi

    if [[ ! ${MERGE} ]]; then
        if [[ ! -e FLAG_${SAMPLE}_MAP ]]; then
            echo - Running bwa mem on ./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz ./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq    
            time bwa mem  ${REF} \
                    ./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz ./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq.gz \
                    | samtools sort -O BAM -@ 16 -o ./03-mapping/${SAMPLE}_sorted.bam \
                    > ./03-mapping/${SAMPLE}_map.log 2>&1

            touch FLAG_${SAMPLE}_MAP 
        fi
    else
        if [[ ! -e FLAG_${SAMPLE}_MAP ]]; then
            echo - Running bwa mem  on ./02-merge/${SAMPLE}_merged.fastq.gz     
            time bwa mem  ${REF}  ./02-merge/${SAMPLE}_merged.fastq.gz \
                    | samtools sort -O BAM -@ 16 -o ./03-mapping/${SAMPLE}_sorted.bam \
                    > ./03-mapping/${SAMPLE}_map.log 2>&1 

            touch FLAG_${SAMPLE}_MAP
        fi
    fi

    #########################################
    # Generate bam index for IGV viewer
    #########################################

    if [[ ! -e FLAG_${SAMPLE}_SORT_BAM ]]; then
        echo - indexing ./03-mapping/${SAMPLE}_sorted.bam     
        samtools index -o ./03-mapping/${SAMPLE}_sorted.bam.bai  ./03-mapping/${SAMPLE}_sorted.bam
        touch FLAG_${SAMPLE}_SORT_BAM
    fi

    #############################################
    # Add RG field 
    #############################################
    if [[ ! -e FLAG_${SAMPLE}_RG  ]]; then
        time java -jar $PICARD_JAR AddOrReplaceReadGroups -I ./03-mapping/${SAMPLE}_sorted.bam \
                -O ./03-mapping/${SAMPLE}_sorted_RG.bam \
                -ID ${SAMPLE} -LB lib1 -PL ILLUMINA -PU unit1 \
                -SM ${SAMPLE} > ./03-mapping/${SAMPLE}_add_group_id.log 2>&1

        samtools index -o ./03-mapping/${SAMPLE}_sorted_RG.bam.bai ./03-mapping/${SAMPLE}_sorted_RG.bam
        touch FLAG_${SAMPLE}_RG
    fi

    #############################################
    # Remove duplicated reads
    #############################################
   
    if [[ ! -e FLAG_${SAMPLE}_REMOVE_DUP  ]]; then
        time gatk --java-options "-Xmx64G -Djava.io.tmpdir=${TMPDIR} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true'" MarkDuplicates \
                -O ./03-mapping/${SAMPLE}_sorted_RG.dedup.bam \
                -I ./03-mapping/${SAMPLE}_sorted_RG.bam \
                -M ./03-mapping/${SAMPLE}_sorted_RG.metrics.txt \
                --REMOVE_DUPLICATES true \
                --TMP_DIR ${TMPDIR}
                
        samtools index -o ./03-mapping/${SAMPLE}_sorted_RG.dedup.bam.bai ./03-mapping/${SAMPLE}_sorted_RG.dedup.bam
        touch FLAG_${SAMPLE}_REMOVE_DUP
    fi
    

    ###########################
    # Run ABRA2 
    ###########################

    if [[ ! -d "04-abra2" ]]; then
        mkdir "04-abra2"
    fi

    if [[ ! -e FLAG_${SAMPLE}_ABRA2 ]]; then
        if [[ ! ${MERGE} ]]; then
            
                time java -Xmx16G -jar abra2-2.23.jar \
                    --in ./03-mapping/${SAMPLE}_sorted_RG.dedup.bam \
                    --out ./04-abra2/${SAMPLE}_realigned.bam \
                    --ref ${REFERENCE} \
                    --threads 16 \
                    --tmpdir ./abra2_tmpdir \
                    --single > abra.log 2>&1

                touch FLAG_${SAMPLE}_ABRA2
            
        else
                time java -Xmx16G -jar abra2-2.23.jar \
                    --in ./03-mapping/${SAMPLE}_sorted_RG.dedup.bam \
                    --out ./04-abra2/${SAMPLE}_realigned.bam \
                    --ref ${REFERENCE} \
                    --threads 16 \
                    --tmpdir ./abra2_tmpdir \
                    > abra.log 2>&1

                samtools index ./04-abra2/${SAMPLE}_realigned.bam

                touch FLAG_${SAMPLE}_ABRA2
        fi
    fi

    #############################################
    # Extract variants and search for INDELs only
    #############################################

    if [[ ! -d "05-variants" ]]; then
        mkdir "05-variants"
    fi


    BAM=./04-abra2/${SAMPLE}_realigned.bam


    # SNP-INDEL calls using GATK
    if [[ ! -e FLAG_${SAMPLE}_GATK_HC  ]]; then
    
        echo; echo "Running gatk HaplotypeCaller -ERC GVCF on ${BAM}"; echo
        time gatk --java-options "-Xmx8G" HaplotypeCaller \
            -R ${REFERENCE} \
            -I ${BAM} \
            -ERC GVCF \
            --stand-call-conf 0  \
            --min-base-quality-score 10 \
            -O ./05-variants/${SAMPLE}.raw.snps.indels.g.vcf \
            -create-output-variant-index true \
            --sample-ploidy 2

    touch FLAG_${SAMPLE}_GATK_HC     
    fi        

done < $PREFIX_FILE

merge_param=()
while read SAMPLE
do
    merge_param+="--variant ./05-variants/${SAMPLE}.raw.snps.indels.g.vcf "
done < $PREFIX_FILE

##################
# Combine gvcf files
##################



if [[ ! -e FLAG_ALL_GATK_COMBINE_GVCFS  ]]; then

    echo; echo "Running GenotypeGVCFs on *.g.vcf files" ; echo
    time gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true'" CombineGVCFs \
        -R ${REFERENCE} \
        -O ./05-variants/ALL_raw_variants_combined.vcf \
        $merge_param
        
        touch FLAG_ALL_GATK_COMBINE_GVCFS
fi

if [[ ! -e FLAG_ALL_GATK_CALL_GVCF  ]]; then

    echo; echo "Running GenotypeGVCFs on *.g.vcf files" ; echo
    time gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true'" GenotypeGVCFs \
        -R ${REFERENCE} \
        --max-alternate-alleles 6 \
        --stand-call-conf 30.0 \
        -O ./05-variants/ALL_raw_variants.vcf \
        -V ./05-variants/ALL_raw_variants_combined.vcf
        
        touch FLAG_ALL_GATK_CALL_GVCF
fi

if [[ ! -e FLAG_ALL_GATK_RAW_VARIANTS ]]; then

    echo; echo "Running SelectVariants on ALL_raw_variants.vcf" ; echo
    time gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR}" SelectVariants \
        -R ${REFERENCE} \
        -V ./05-variants/ALL_raw_variants.vcf \
        --select-type-to-include SNP \
        --select-type-to-include INDEL \
        -O ./05-variants/ALL_raw_snps_indels.vcf

    touch FLAG_ALL_GATK_RAW_VARIANTS
fi


if [[ ! -e FLAG_ALL_GATK_FILTERED_VARIANTS ]]; then
    echo; echo "Running VariantFiltration on ALL_raw_snps.vcf" ; echo
    time gatk --java-options "-Xmx8G -Djava.io.tmpdir=${TMPDIR}" VariantFiltration \
        -R ${REFERENCE} \
        -V ./05-variants/ALL_raw_snps_indels.vcf \
        --filter-name  "filter_QD" --filter-expression "QD < 2.0" \
        --filter-name  "filter_FS" --filter-expression "FS > 60.0" \
        --filter-name  "filter_MQ" --filter-expression "MQ < 40.0" \
        --filter-name  "filter_MQRankSum" --filter-expression "MQRankSum < -12.5" \
        --filter-name  "ftr_ReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
        -O ./05-variants/ALL_filtered_snps_indels.vcf

    touch FLAG_ALL_GATK_FILTERED_VARIANTS
fi


#################
# Annotate InDels
#################

if [[ ! -e FLAG_ALL_ANNOTATE_VARIANTS ]]; then
    echo; echo Annotating ./05-variants/ALL_filtered_snps_indels.vcf; echo
    java -Xmx30g -jar $SNPEFF_JAR  -c ./snpEff.config -v -lof -motif -hgvs -nextProt GRCm38.99 ./05-variants/ALL_filtered_snps_indels.vcf > ./05-variants/ALL_filtered_snps_indels.snpEff.vcf    

    touch FLAG_ALL_ANNOTATE_VARIANTS
fi

