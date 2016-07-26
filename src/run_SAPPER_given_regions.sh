#!/bin/bash

# This is a template pipeline script to call SNVs from any given
# target regions from ChIP and ctonrol BAM files.

# Please read the descriptions and modify the following values for
# your own datasets.

# NOTE: ALL INPUT FILES SHOULD BE AT THE CURRENT WORKING DIRECTORY

# Name of the factor or the experiment
FACTORNAME="MyFactor";

# Path to the IP BAM file
CHIPBAM="IP.bam";

# Path to the control BAM file
CTRLBAM="CTRL.bam";

# Path to the region file. Must be BED-style file containing 1st
# column as chromosome name, 2nd column as start and 3rd as end
# position
BEDFILE="Target.bed";

# Filtering criteria for SAPPER

# Minimum depth for SAPPER
MINDEPTH=20

# Minimum Genotype Quality Score cutoff for heterozygous loci
MINHETEROGQ=100

# Minimum Genotype Quality Score cutoff for homozygous loci
MINHOMOGQ=10

###################DO NOT EDIT CODES BELOW################

# check settings

# step 1 filter
samtools view -q 30 -F 4 -F 256 -F 2048 -b ${CHIPBAM} -o ${CHIPBAM/.bam/_clean.bam} &
samtools view -q 30 -F 4 -F 256 -F 2048 -b ${CTRLBAM} -o ${CTRLBAM/.bam/_clean.bam} &

wait;

# step 2 sort
samtools sort  ${CHIPBAM/.bam/_clean.bam} ${CHIPBAM/.bam/_clean_sorted} &
samtools sort  ${CTRLBAM/.bam/_clean.bam} ${CTRLBAM/.bam/_clean_sorted} &

wait;

# step 3 sort bed file
sort -k1,1 -k2,2n ${BEDFILE} > ${FACTORNAME}.sorted.bed &

wait;

# step 4 extract reads in target regions towards each side for 150bps then merge:

awk -v OFS="\t" '{print $1,(($2-150)>0?($2-150):0),$3+150}' ${FACTORNAME}.sorted.bed | bedtools merge -i - > ${FACTORNAME}_extended.bed


samtools view -b ${CHIPBAM/.bam/_clean_sorted.bam}  -L ${FACTORNAME}_extended.bed -o ${FACTORNAME}_CHIP_regions.bam &
samtools view -b ${CTRLBAM/.bam/_clean_sorted.bam}  -L ${FACTORNAME}_extended.bed -o ${FACTORNAME}_CTRL_regions.bam &

wait;

# step 5 run SAPPER
SAPPER call -b ${FACTORNAME}_extended.bed -t ${FACTORNAME}_CHIP_regions.bam -c ${FACTORNAME}_CTRL_regions.bam -o ${FACTORNAME}_SAPPER.vcf &

wait;

# step 6 run SAPPER_filter

SAPPER filter -i ${FACTORNAME}_SAPPER.vcf -d ${MINDEPTH} -t homo -q ${MINHOMOGQ} -o ${FACTORNAME}_SAPPER_homo.vcf
SAPPER filter -i ${FACTORNAME}_SAPPER.vcf -d ${MINDEPTH} -t hetero -q ${MINHETEROGQ} -o ${FACTORNAME}_SAPPER_hetero.vcf

echo "Finished! Please check:"
echo " " ${FACTORNAME}_SAPPER.vcf "for all the predicted SNVs."
echo " " ${FACTORNAME}_SAPPER_hetero.vcf "for the filtered heterozygous SNVs above cutoff."
echo " " ${FACTORNAME}_SAPPER_homo.vcf "for the filtered homozygous SNVs above cutoff."

echo "Please use SAPPER_filter script to further filter the results."
