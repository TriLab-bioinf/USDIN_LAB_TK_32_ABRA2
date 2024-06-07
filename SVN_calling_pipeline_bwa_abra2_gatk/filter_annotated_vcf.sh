#!/usr/bin/bash

VCF=$1
FILTER_FILE=$2

FILTER_OUT=()
while read filter; do
    FILTER_OUT+="${filter}\|"
    
done < $FILTER_FILE

# Removing last 2 characters "\|"
FILTER_OUT=${FILTER_OUT::-2}

echo "#Filter_out = <<${FILTER_OUT}>>";echo
perl -e 'print "#CHROM\tPOS\tTYPE\tGENE\tEFFECT\tGAP_LENGTH\tREAD_DEPTH\n"'

#grep -v 'intergenic\|intron_variant\|downstream_gene_variant\|upstream_gene_variant' ${VCF} | \
grep -vw "${FILTER_OUT}" ${VCF} | \
perl -lane 'chomp;
            $gap = (length($F[3]) - length($F[4]));
            @ann = split(/\|/, $F[7]);
            $idv = $1 if $F[7] =~ m/IDV=(\d+)/;
            @sample1 = split(":",$F[9]);
            @sample2 = split(":",$F[10]);
            @sample3 = split(":",$F[11]);

            $sample1[0] =~ s:[/|]:_:;
            $sample2[0] =~ s:[/|]:_:;
            $sample3[0] =~ s:[/|]:_:;

            #print "$sample1[0] ; $sample2[0] ; $sample3[0]";
            if ( ($sample1[0] ne $sample2[0] || $sample1[0] ne $sample3[0] ) && $sample1[0] ne "._."){
            
                print join("\t", ($F[0], $F[1], $F[3], $F[4], $ann[1], $ann[3], $ann[7], $ann[10], $gap, $idv, $F[8], $F[9], $F[10], $F[11]));
        
            }
        '
