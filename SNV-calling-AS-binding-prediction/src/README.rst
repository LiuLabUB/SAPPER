Install
=======

Simply go to SNVAS source code dir then compile SNVAS and the modified
fermi (https://github.com/lh3/fermi/)::

 $ make

You will have executable binary files 'fermi', 'SNVAS' and
'SNVAS_filter' in the 'bin' subdirectory. Now copy/move them to one of
your PATH such as /usr/loca/bin::

 $ mv fermi /usr/local/bin
 $ mv SNVAS /usr/local/bin
 $ mv SNVAS_filter /usr/local/bin

Usage
=====

Before running SNVAS/Preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. SNVAS only accepts base quality scores with Illumina 1.8+
   (Phred+33) coding. For Illumina 1.3+ or 1.5+ (Phred+64) coding, the
   fastq file should be changed to Phred+33 type. We recommand using the
   tool "seqtk" (https://github.com/lh3/seqtk).

2. Map fastq file to reference genome, and get the BAM file. For
   paired-end data, we recommend 'bwa mem', and for single-end, 'bwa
   aln'.

3. Remove reads flagged as low mapping quality (MAPQ), unmapped,
   secondary alignment or supplementary alignment in the BAM file. Check
   both ChIP-seq and control SAM/BAM files. Such as::

      $ samtools view -q 30 -F 4 -F 256 -F 2048 -bS sample.sam -o sample_filter.bam

4. Sort the BAM file after filtering accroding to coordinate, using
   samtools or Picard::

      $ samtools sort  sample_filter.bam  sample_filter_sorted

5. Peak calling. We recommand using the software "macs2"
   (https://github.com/taoliu/MACS).

   Exampel for paired-end ChIP-seq::

      $ macs2 callpeak -f BAMPE -t CHIP_filtered_sorted.bam -c Ctrl_filtered_sorted.bam -n MyFactor -g hs


   And for single-end ChIP-seq::

      $ macs2 callpeak -f BAM -t CHIP_filtered_sorted.bam -c Ctrl_filtered_sorted.bam -n MyFactor -g hs

   Then the peak region file 'MyFactor_peaks.narrowPeak' should be used
   in downstream analysis. Then the peak regions should be sorted
   according to coordinates::

      $ sort -k1,1 -k2,2n MyFactor_peaks.narrowPeak > MyFactor_peaks.sorted.bed

6. Extract reads in selected peak region, and generate the subset BAM
   files from both ChIP-seq and control dataset. Such as::

      $ samtools view -b sample_filter_sorted.bam -L sample_peaks_sorted.bed -o sample_peaks_sorted.bam

Finally, there are 3 files which should be prepared before running
SNVAS: a bed file for peak regions; a BAM file of ChIP-seq reads
mapped to the peak regions; and a BAM file of control reads mapped 
to the peak region. All of them are sorted by genome coordinate.

Running SNVAS
~~~~~~~~~~~~~

1. To get a listing of all parameters, run ``SNVAS -h``.

2. For paired-end data, you can run::

     $ SNVAS sample_peaks_sorted.bed sample_peaks_sorted.bam control_peaks_sorted.bam PE sample.vcf

   PE is the parameter shows the data is paired-end. sample.vcf is the
   output vcf file

3. For single-end data, you should change "PE" to "SE".

Batch script
~~~~~~~~~~~~

We provided a shell script ``run_SNVAS.sh`` to go through the above
steps in a pipeline. Please open and edit the ``run_SNVAS.sh`` file,
and put it in the working directory where your BAM files are
located. Please check the description in the shell script.

Interpret Results
=================

Format of the output
~~~~~~~~~~~~~~~~~~~~

The output vcf file follows VCF v4.1 format
(http://samtools.github.io/hts-specs/VCFv4.1.pdf); the detail
information of each term is defined in the header of the vcf file::

 ##fileformat=VCFv4.1
 ##fileDate=2015817
 ##source=SNVAS_V0.1
 ##Program Args: sample_peaks_sorted.bed sample_peaks_sorted.bam control_peaks_sorted.bam PE sample.vcf
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
   a heterozygous SNV "MinBIC_model:homo", 2) heterozygous SNV with
   allele specific binding "MinBIC_model:heter_AS", or 3) heterozygous
   SNV without allele specific binding from our model
   "MinBIC_model:heter_noAS".

4. We use genotype quality score to measure the reliability of the
   predicted SNVs. For the homozygous SNV, see the term "GQ_homo"; for
   the allele-specifically bound heterozygous SNV, see the term
   "GQ_heter_AS"; for the non allele-specifically bound heterozygous SNV,
   see the term "GQ_heter_noAS". Higher the genotype quality score,
   more reliable the prediction is. 

Note, there is no cutoff applied in the VCF file. The only rule is the
BIC, so that the reported genotype/allele-specific status has the
smallest BIC among all the other models. We provide downstream
analysis tool 'SNVAS_filter' to further filter the results in VCF
files.

Filtering results using SNVAS_filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We provided a postprocessing tool ``SNVAS_filter`` to further filter
the output VCF file. It can be used to get a list of 1) homozygous
SNVs; 2) heterozygous SNVs; 3) heterozygous SNVs with non-allele
specific binding; 4) heterozygous SNVs with allele-specific binding:

1. To get a listing of all parameters, run ``SNVAS_filter -h``.

2. To get homozygous SNVs with quality score >=cutoff (integer), you
   can run::

      $ SNVAS_filter sample.vcf homo cutoff sample_homo_afterfilter.vcf

3. To get all heterozygous SNVs with quality score >=cutoff (integer),
   you can run::

      $ SNVAS_filter sample.vcf hete cutoff sample_hete_afterfilter.vcf

4. To get allele-specific heterozygous SNVs with quality score
   >=cutoff (integer), you can run::

      $ SNVAS_filter sample.vcf heter_AS cutoff sample_heterAS_afterfilter.vcf

5. To get non allele-specific heterozygous SNV with quality score
   >=cutoff (integer), you can run::

      $ SNVAS_filter sample.vcf heter_noAS cutoff sample_heterNonAS_afterfilter.vcf


Release Notes
=============
Release 0.1 (2015-08-14)
This is the first public release of SNVAS.
