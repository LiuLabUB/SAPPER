#!/bin/bash

# Please read the descriptions and modify the following values for
# your own datasets.

# NOTE: ALL INPUT FILES SHOULD BE AT THE CURRENT WORKING DIRECTORY

# Name of the factor or the experiment
FACTORNAME="MyFactor";

# Path to the IP BAM file
CHIPBAM="IP.bam";

# Path to the control BAM file
CTRLBAM="CTRL.bam";

# For PE data, use 'PE', or for SE data, use 'SE'.
MODE="PE";

# Extra parameters for MACS2. Note, do not change '--broad' setting
# here. Do such in the next setting for MACS2MODE. And do not set '-f
# BAMPE' here if you are dealing with PE data, it will be set
# according to your previous MODE setting. If you plan to use "PE"
# MODE, you have to set '--broad-cutoff' together with '-q' or '-p'.
MACS2EXTPARAM="-g hs -q 0.05 --fe-cutoff 5";
#MACS2EXTPARAM="-g hs --broad-cutoff 0.1 -p 0.05 --fe-cutoff 5"; # this is for broad peak calling

# For broad region calling, use 'broad'; for narrow peak calling, use
# 'narrow'. 
MACS2MODE="narrow";

###################DO NOT EDIT CODES BELOW################

# check settings

if [ ${MACS2MODE} = "narrow" ];then
    MACS2OUTPUTFILE="${FACTORNAME}_peaks.narrowPeak";
elif [ ${MACS2MODE} = "broad" ];then
    MACS2OUTPUTFILE="${FACTORNAME}_peaks.broadPeak";
    MACS2EXTPARAM="${MACS2EXTPARAM} --broad";
else
    echo "Wrong MACS2MODE setting! Please choose from 'narrow' and 'broad'!";
    exit;
fi

if [ ${MODE} = "PE" ];then
    MACS2EXTPARAM="${MACS2EXTPARAM} -f BAMPE";
elif [ ${MODE} = "SE" ];then
    MACS2EXTPARAM="${MACS2EXTPARAM} -f BAM";
else
    echo "Wrong MODE setting! Please choose from 'PE' and 'SE'!";
    exit;
fi  

# step 1 filter
samtools view -q 30 -F 4 -F 256 -F 2048 -b ${CHIPBAM} -o ${CHIPBAM/.bam/_clean.bam} &
samtools view -q 30 -F 4 -F 256 -F 2048 -b ${CTRLBAM} -o ${CTRLBAM/.bam/_clean.bam} &

wait;

# step 2 sort
samtools sort  ${CHIPBAM/.bam/_clean.bam} ${CHIPBAM/.bam/_clean_sorted} &
samtools sort  ${CTRLBAM/.bam/_clean.bam} ${CTRLBAM/.bam/_clean_sorted} &

wait;

# step 3 peak calling
macs2 callpeak -t ${CHIPBAM/.bam/_clean_sorted.bam}  -c ${CTRLBAM/.bam/_clean_sorted.bam} -n ${FACTORNAME} ${MACS2EXTPARAM} &

wait;

# step 4 sort peak file
sort -k1,1 -k2,2n ${MACS2OUTPUTFILE} > ${FACTORNAME}_peaks.sorted.bed &

wait;

# step 5 extract reads in peak regions towards each side for 150bps then merge:

awk -v OFS="\t" '{print $1,(($2-150)>0?($2-150):0),$3+150}' ${FACTORNAME}_peaks.sorted.bed | bedtools merge -i - > ${FACTORNAME}_extended_peaks.bed


samtools view -b ${CHIPBAM/.bam/_clean_sorted.bam}  -L ${FACTORNAME}_extended_peaks.bed -o ${FACTORNAME}_CHIP_peaks.bam &
samtools view -b ${CTRLBAM/.bam/_clean_sorted.bam}  -L ${FACTORNAME}_extended_peaks.bed -o ${FACTORNAME}_CTRL_peaks.bam &

wait;

# step 6 run SNVAS
SNVAS ${FACTORNAME}_extended_peaks.bed ${FACTORNAME}_CHIP_peaks.bam ${FACTORNAME}_CTRL_peaks.bam ${FACTORNAME}_SNVAS.vcf &

wait;

echo "Finished. Please check the file" ${FACTORNAME}_SNVAS.vcf
echo "Please use SNVAS_filter script to further filter the results."
