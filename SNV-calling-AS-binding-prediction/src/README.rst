Install
=======

Compile 'fermi' https://github.com/lh3/fermi/
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

```
$ cd fermi
$ make
$ cd ..
```

You will have executable fermi binary file under 'fermi'
directory. Now copy/move it to one of your PATH such as
/usr/local/bin.

```
$ mv fermi /usr/local/bin
```


Then Compile SNVAS
~~~~~~~~~~~~~~~~~~

```
$ make
```

You will have executable binary files 'SNVAS' under the current src
directory. Now copy/move it to one of your PATH such as /usr/loca/bin.

```
$ mv SNVAS /usr/local/bin
```

Usage
=====

Before running SNVAS/Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. SNVAS only accepts base quality scores with Illumina 1.8+ (Phred+33) coding. For Illumina 1.3+ or 1.5+ (Phred+64) coding, the fastq file should be changed to Phred+33 type. We recommand using the tool "seqtk" (https://github.com/lh3/seqtk).

2. Map fastq file to reference genome, and get the BAM file. For
paired-end data, we recommend 'bwa mem', and for single-end, 'bwa aln'.

3. Remove reads flagged as low mapping quality (MAPQ), unmapped,
secondary alignment or supplementary alignment in the BAM file. Check
both ChIP-seq and control SAM/BAM files.

e.g. 

```
$ samtools view -q 30 -F 4 -F 256 -F 2048 -bS sample.sam -o sample_filter.bam
```

4. Sort the BAM file after filtering accroding to coordinate, using
samtools or Picard.

```
$ samtools sort  sample_filter.bam  sample_filter_sorted
```

5. Peak calling. We recommand using the software "macs2"
(https://github.com/taoliu/MACS).

Exampel for paired-end ChIP-seq:

```
$ macs2 callpeak -f BAMPE -t CHIP_filtered_sorted.bam -c Ctrl_filtered_sorted.bam -n MyFactor -g hs
```

And for single-end ChIP-seq:

```
$ macs2 callpeak -f BAM -t CHIP_filtered_sorted.bam -c Ctrl_filtered_sorted.bam -n MyFactor -g hs
```

Then the peak region file 'MyFactor_peaks.narrowPeak' should be used
in downstream analysis. Then the peak regions should be sorted
according to coordinates.

```
sort -k1,1 -k2,2n MyFactor_peaks.narrowPeak > MyFactor_peaks.sorted.bed
```


6. Extract reads in selected peak region, and generate the subset BAM files
from both ChIP-seq and control dataset.

e.g.

```
samtools view -b sample_filter_sorted.bam -L sample_peaks_sorted.bed -o sample_peaks_sorted.bam
```

Finally, there are 3 files which should be prepared before running SNVAS

1. peak region bed file (sorted by coordinate)

2. bam file of ChIP-seq dataset in the peak region (sorted by coordinate) 

3. bam file of control dataset in the peak region (sorted by coordinate)

Running SNVAS
~~~~~~~~~~~~~

1. To get a listing of all parameters, run ```SNVAS -h```.

2. For paired-end data, you can run:

```
$ SNVAS sample_peaks_sorted.bed sample_peaks_sorted.bam control_peaks_sorted.bam PE sample.vcf
```

PE is the parameter shows the data is paired-end. sample.vcf is the output vcf file

3. For single-end data, you should change "PE" to "SE".

Interpret Results
=================
The output vcf file is based on VCF v4.1 format; the detail information of each term is defined in the header of the vcf file.

Note:
1 Column1(CHROM) and column2(POS) define the position (1-based) of the SNV.

2 The last column stands for the SNV type. "0|1" or "1|2" stands for heterozygous SNV, and "1|1" stands for homozygous SNV. By now, this tool is only for single ChIP-seq data, so there is no "0|2", "2|2" or other type.

3 The term "MinBIC_model" defines the heterozygous SNV is allele-specific binding or not. "MinBIC_model:heter_AS" means it is allele-specific, while "MinBIC_model:heter_noAS" means it is not allele-specific; "MinBIC_model:homo" means it is homozygous SNV.

4 We use genotype quality score to measure the reliability of predicted SNV. For homozygous SNV, see the term "GQ_homo"; for allele-specific heterozygous SNV, see the term "GQ_heter_AS"; for no allele-specific heterozygous SNV, see the term "GQ_heter_noAS". The genotype quality score is higher, then the SNV is more reliable.

Release Notes
=============
Release 0.1 (2015-08-12)
This is the first public release of SNVAS.
