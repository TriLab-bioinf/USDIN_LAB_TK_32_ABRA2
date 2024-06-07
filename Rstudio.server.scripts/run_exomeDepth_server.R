#!/usr/bin/env Rscript
library("optparse")
library("dplyr")

# Called from the following commands:
# Rscript ./run_exomeDepth.R -i ${SAMPLE}_realigned_nochrM.bam -c ${CONTROL_BAMS} -o 07-exomeDepth -r ${REFERENCE}

# Parse options 
option_list = list(
    make_option(c("-i", "--input_bam_file"), type="character", default=NULL, 
              help="input bam file", metavar="character"),
    make_option(c("-c", "--control_bam_file"), type="character", default=NULL, 
              help="file containing names of control bam files", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default="ExomeDepth_DIR", 
              help="ExomeDepth output directory [default= %default]", metavar="character"),
    make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="genome reference fasta file", metavar="character")
    )
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load libraries
suppressMessages(library("TxDb.Mmusculus.UCSC.mm10.knownGene"))
suppressMessages(library("plyranges"))
suppressMessages(library("ExomeDepth"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library("AnnotationHub"))
suppressMessages(library("AnnotationDbi"))

# Load main exomedepth function
source(file = "aux1-exomedepth.R")

####################################################################
# Set inputs
####################################################################

# Upload control bam files (to contrast with sample bam)
#control_bams <- read.delim(file = opt$control_bam_file, header = FALSE)
control_bams <- read.delim(file = "./control_bams.txt", header = FALSE)

#input_bam_file  <-  opt$input_bam_file # "LB2BAV7596_realigned_nochrM.bam"
#reference <- opt$reference #"genome.fa"

case_files <- c("Slx4-A24_realigned_nochrM.bam","Slx4-A25_realigned_nochrM.bam","LB2BAV7596_realigned_nochrM.bam")

# Choose one input sample at a time by changing the index of case_files:
input_bam_file  <-   case_files[3]

reference <- "genome.fa"
####################################################################

# Load annotation data from UCSC
txdb <- txdbmaker::makeTxDbFromGFF("/fdb/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf") #TxDb.Mmusculus.UCSC.mm10.knownGene
mykeys <- keys(txdb, keytype = "GENEID")


# Create datafram with columns in the correct order (chr, start, end, others....)
exons.mm10 <- AnnotationDbi::select(txdb,
                                    keys = mykeys,
                                    keytype="GENEID",
                                    columns = c("EXONCHROM", "EXONSTART", "EXONEND","GENEID", "EXONID")
                                    )[,c(3,4,5,1,2)]

# exons.mm10 <- AnnotationDbi::select(txdb,
#        keys = mykeys,
#        keytype="GENEID",
#        columns = c("EXONCHROM", "EXONSTART", "EXONEND","GENEID", "EXONID")
#       )[,c(2,3,4,1,5)]

# Remove non-chromosomal contigs and mitochondrial chromosome. chrX and chrY should be removed as well if the controls have a different gender to the actual sample
keep <- grep(pattern = '^(chr[0-9]+|chrX|chrY)$', exons.mm10$EXONCHROM)
exons.mm10 <- exons.mm10[keep,]

# Generate GR object with annotation (used later...)
exons.mm10.GRanges <- GenomicRanges::GRanges(seqnames = exons.mm10$EXONCHROM,
                                                IRanges::IRanges(
                                                 start=exons.mm10$EXONSTART,
                                                 end=exons.mm10$EXONEND
                                                 ),
                                                names = exons.mm10$GENEID,
                                                exon_id = exons.mm10$EXONID
                                            )



myBams <- c(paste0('./04-abra2/',input_bam_file), control_bams[,1])

ExomeCount <- getBamCounts(bed.frame=exons.mm10,
                           bam.files=myBams,
                           include.chr=F,
                           referenceFasta=reference
                           )



### Identify the reference set of bam files (note the name change due to importing into R)
my.ref.samples <- c(control_bams[,1])

# Remove bam extension from input file name
sample_name <- gsub(".bam","",input_bam_file)

#outdir <- opt$outdir #"07-exomeDepth"
outdir <- "07-exomeDepth"
dir.create(path = outdir, showWarnings = FALSE)

# Run ExomeDepth

print('Running exomeDepth ...')

exomeDepth(my.test.num = input_bam_file,
        outdir = outdir,
        ExomeCount = ExomeCount,
        my.ref.samples = my.ref.samples,
        exons.mm10.GRanges = exons.mm10.GRanges,
        sample=sample_name
        #genes.mm10.GRanges=genes.mm10.GRanges
        ) 

print(paste("Analysis complete for", input_bam_file))
