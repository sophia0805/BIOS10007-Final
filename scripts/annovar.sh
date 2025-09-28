#!/bin/bash

invcf=$1
prefix=$2

annovarpath=../software/annovar
humandb=../software/annovar/humandb

##########################################
# STEP 1 - Convert vcf to annovar format #
##########################################
$annovarpath/convert2annovar.pl -format vcf4 $invcf -outfile ${prefix}.annovar

##################################
# STEP 2 - Annotate the variants #
##################################

# Annotate with refgene info
#$annovarpath/annotate_variation.pl --geneanno --buildver hg38 -dbtype refGene ${prefix}.avinput $humandb

# Annotate with GWAS info
#$annovarpath/annotate_variation.pl -regionanno --buildver hg38 -dbtype clinvar_20221231 ${prefix}.avinput $humandb

# Annotate with dbSNP info
#$annovarpath/annotate_variation.pl -regionanno --buildver hg38 -dbtype avsnp150 ${prefix}.avinput $humandb

# Annotate with various functional scores
#$annovarpath/annotate_variation.pl -regionanno --buildver hg38 -dbtype dbnsfp42c ${prefix}.avinput $humandb

# table_annovar does all four commands at once and outputs a file that is easily usable in Excel
$annovarpath/table_annovar.pl ${prefix}.annovar $humandb --buildver hg38 --protocol refGene,avsnp150,clinvar_20221231,dbnsfp42c --operation g,f,f,f
