Install
=======

Simply go to SAPPER source code dir then compile SAPPER and the modified
fermi (https://github.com/lh3/fermi/)::

 $ make

You will have executable binary files 'SAPPER' in the current directory. Now
copy/move them to one of your PATH such as /usr/loca/bin::

 $ mv SAPPER /usr/local/bin

Usage
=====

Before running SAPPER/Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. SAPPER only accepts base quality scores with Illumina 1.8+
   (Phred+33) coding. For Illumina 1.3+ or 1.5+ (Phred+64) coding, the
   fastq file should be changed to Phred+33 type. We recommand using the
   tool "seqtk" (https://github.com/lh3/seqtk).

2. Map fastq files to reference genome, and get the BAM files. For
   paired-end data, we recommend 'bwa mem', and for single-end, 'bwa
   aln'.

3. Remove reads flagged as low mapping quality (MAPQ), unmapped,
   secondary alignment or supplementary alignment in the BAM file. Check
   both ChIP-seq and control SAM/BAM files. Such as::

      $ samtools view -q 30 -F 4 -F 256 -F 2048 -bS sample.sam -o sample_filter.bam

4. Sort the BAM file after filtering accroding to coordinates, using
   samtools or Picard::

      $ samtools sort  sample_filter.bam  sample_filter_sorted

5. Peak calling. We recommand using MACS2 (https://github.com/taoliu/MACS).

   Example for paired-end ChIP-seq::

      $ macs2 callpeak -f BAMPE -t CHIP_filtered_sorted.bam -c Ctrl_filtered_sorted.bam -n MyFactor -g hs

   And for single-end ChIP-seq::

      $ macs2 callpeak -f BAM -t CHIP_filtered_sorted.bam -c Ctrl_filtered_sorted.bam -n MyFactor -g hs

   Optionally, you can ask MACS2 to apply a fold-enrichment cutoff
   (--fe-cutoff) so that only high quality peaks are kept in
   downstream analysis.

   The peak region file 'MyFactor_peaks.narrowPeak' should then be
   sorted according to coordinates::

      $ sort -k1,1 -k2,2n MyFactor_peaks.narrowPeak > MyFactor_peaks.sorted.bed

6. Extract reads in selected peak region with extension to each side
   (e.g. 300bps), and generate the subset BAM files from both ChIP-seq
   and control dataset. Such as::

      $ awk -v OFS="\t" '{print $1,(($2-150)>0?($2-150):0),$3+150}' MyFactor_peaks.sorted.bed | bedtools merge -i - > MyFactor_extended_peaks.bed
      $ samtools view -b sample_filter_sorted.bam -L MyFactor_extended_peaks.bed -o reads_in_extended_peaks.bam

Finally, there are 3 files which should be prepared before running
SAPPER: a bed file for peak regions; a BAM file of ChIP-seq reads
mapped to the peak regions; and a BAM file of control reads mapped 
to the peak region. All of them are sorted by genome coordinate.

Running SAPPER
~~~~~~~~~~~~~

1. Simply run::

     $ SAPPER call -b sample_peaks_sorted.bed -t sample_peaks_sorted.bam -c control_peaks_sorted.bam -o sample.vcf

And the sample.vcf is the output. SAPPER tool outputs all possible
SNVs. We recommend to use the filter tool (see below) to get a
reliable list of SNVs above certain cutoff values.

Batch script
~~~~~~~~~~~~

We provided a shell script ``run_SAPPER.sh`` to go through the above
steps, and additional filtering in a pipeline. Please open and edit
the ``run_SAPPER.sh`` file, and put it in the working directory where
your BAM files are located. Please check the description in the shell
script.

Interpret Results
=================

Format of the output
~~~~~~~~~~~~~~~~~~~~

The output vcf file follows VCF v4.1 format
(http://samtools.github.io/hts-specs/VCFv4.1.pdf); the detail
information of each term is defined in the header of the vcf file::

 ##fileformat=VCFv4.1
 ##fileDate=2015817
 ##source=SAPPER_V0.1
 ##Program Args="sample_peaks_sorted.bed sample_peaks_sorted.bam control_peaks_sorted.bam PE sample.vcf"
 ##INFO=<ID=MinBIC_model,Number=.,Type=String,Description="Model with minimum BIC value">
 ##INFO=<ID=DP_ChIP,Number=1,Type=Integer,Description="Approximate read depth in ChIP-seq data; some reads may have been filtered">
 ##INFO=<ID=DP_input,Number=1,Type=Integer,Description="Approximate read depth in input data; some reads may have been filtered">
 ##INFO=<ID=fermiNTs,Number=.,Type=String,Description="Nucleotides from the genotype information of fermi assembly result">
 ##INFO=<ID=top1,Number=.,Type=String,Description="Read depth of top1 nucleotide in ChIP-seq data; some reads may have been filtered">
 ##INFO=<ID=top2,Number=.,Type=String,Description="Read depth of top2 nucleotide in ChIP-seq data; some reads may have been filtered">
 ##INFO=<ID=top1input,Number=.,Type=String,Description="Read depth of top1 nucleotide in input data; some reads may have been filtered">
 ##INFO=<ID=top2input,Number=.,Type=String,Description="Read depth of top2 nucleotide in input data; some reads may have been filtered">
 ##INFO=<ID=top1raw,Number=.,Type=Integer,Description="Read depth of top1 nucleotide in raw ChIP-seq data">
 ##INFO=<ID=top2raw,Number=.,Type=Integer,Description="Read depth of top2 nucleotide in raw ChIP-seq data">
 ##INFO=<ID=top1inputraw,Number=.,Type=Integer,Description="Read depth of top1 nucleotide in raw input data">
 ##INFO=<ID=top2inputraw,Number=.,Type=Integer,Description="Read depth of top1 nucleotide in raw input data">
 ##INFO=<ID=lnL_homo_major,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with major allele model">
 ##INFO=<ID=lnL_homo_minor,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with minor allele model">
 ##INFO=<ID=lnL_heter_noAS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with no allele-specific model">
 ##INFO=<ID=lnL_heter_AS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with allele-specific model">
 ##INFO=<ID=BIC_homo_major,Number=1,Type=Float,Description="BIC value of homozygous with major allele model">
 ##INFO=<ID=BIC_homo_minor,Number=1,Type=Float,Description="BIC value of homozygous with minor allele model">
 ##INFO=<ID=BIC_heter_noAS,Number=1,Type=Float,Description="BIC value of heterozygous with no allele-specific model">
 ##INFO=<ID=BIC_heter_AS,Number=1,Type=Float,Description="BIC value of heterozygous with allele-specific model">
 ##INFO=<ID=GQ_homo,Number=1,Type=Float,Description="Genotype quality of homozygous with major allele model">
 ##INFO=<ID=GQ_heter_noAS,Number=1,Type=Float,Description="Genotype quality of heterozygous with no allele-specific model">
 ##INFO=<ID=GQ_heter_AS,Number=1,Type=Float,Description="Genotype quality of heterozygous with allele-specific model">
 ##INFO=<ID=GQ_heter_ASsig,Number=1,Type=Float,Description="Genotype quality of allele-specific significance compared with no allele-specific model">
 ##INFO=<ID=Allele_ratio_heter_AS,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
 ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
 #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE

Important information in the file:

1. The column 1 (CHROM) and column 2 (POS) define the position
   (1-based) of the variation.

2. The last column (SAMPLE) shows the SNV type. "0|1" or "1|2" stands
   for heterozygous SNV, and "1|1" stands for homozygous SNV. By now,
   this tool is only for single ChIP-seq data, so there is no "0|2",
   "2|2" or other type.

3. The term "MinBIC_model" defines the best model (with the smallest
   BIC -- Bayesian Information Criteria) that our method chooses from 1)
   a heterozygous SNV "MinBIC_model=homo", 2) heterozygous SNV with
   allele specific binding "MinBIC_model=heter_AS", or 3) heterozygous
   SNV without allele specific binding from our model
   "MinBIC_model=heter_noAS".

4. We use genotype quality score to measure the reliability of the
   predicted SNVs. For the homozygous SNV, see the term "GQ_homo"; for
   the allele-specifically bound heterozygous SNV, see the term
   "GQ_heter_AS"; for the non allele-specifically bound heterozygous SNV,
   see the term "GQ_heter_noAS". Higher the genotype quality score,
   more reliable the prediction is. 

Note, there is no cutoff applied in the VCF file. The only rule is the
BIC, so that the reported genotype/allele-specific status has the
smallest BIC among all the other models. We provide downstream
analysis tool 'SAPPER_filter' to further filter the results in VCF
files.

Filtering results using SAPPER_filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We provided a postprocessing tool ``SAPPER_filter`` to further filter
the output VCF file. It can be used to get a list of 1) homozygous
SNVs; 2) heterozygous SNVs; 3) heterozygous SNVs with non-allele
specific binding; 4) heterozygous SNVs with allele-specific binding:

1. To get homozygous SNVs::

      $ SAPPER filter -i sample.vcf -d MINDEPTH -t homo -q MINCUTOFF -o sample_homo_afterfilter.vcf

2. To get heterozygous SNVs::

      $ SAPPER filter -i sample.vcf -d MINDEPTH -t hetero -q MINCUTOFF -o sample_hete_afterfilter.vcf

3. To get allele-specific heterozygous SNVs::

      $ SAPPER filter -i sample.vcf -d MINDEPTH -t heter_AS -q MINCUTOFF -o sample_heterAS_afterfilter.vcf

4. To get non allele-specific heterozygous SNV::

      $ SAPPER filter -i sample.vcf -d MINDEPTH -t heter_noAS -q MINCUTOFF -o sample_heterNonAS_afterfilter.vcf

The selection of minimum depth and minimum genotype quality score
cutoffs is arbitrary. We recommand minimum depth of 20, and minimum GQ
100 for heterozygous SNVs and 10 for homozygous SNVs.


Release Notes
=============
Release 0.2 (2016-06-14)
This is the first public release of SAPPER.
