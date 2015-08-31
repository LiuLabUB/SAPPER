Running example
=================

Paired-end ChIP-seq dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Simply run::

     cd PE_demo
     SNVAS call -b PEsample_peaks_sorted.bed -t PEsample_peaks_sorted.bam -c PEcontrol_peaks_sorted.bam -o PEsample.vcf

2. Filtering results::

     SNVAS filter -i PEsample.vcf -o hetero.vcf -t hetero -d 20 -q 100
     SNVAS filter -i PEsample.vcf -o hetero_AS.vcf -t hetero_AS -d 20 -q 100
     SNVAS filter -i PEsample.vcf -o hetero_nonAS.vcf -t hetero_nonAS -d 20 -q 100
     SNVAS filter -i PEsample.vcf -o homo.vcf -t homo -d 20 -q 10

3. Summarize predicted snvs number per kb and ts/tv ratio for different GQ cutoff::

     SNVAS stat -i PEsample.vcf -b PEsample_peaks_sorted.bed -t hetero -o hetero_stat.txt
     SNVAS stat -i PEsample.vcf -b PEsample_peaks_sorted.bed -t hetero_AS -o hetero_AS_stat.txt
     SNVAS stat -i PEsample.vcf -b PEsample_peaks_sorted.bed -t hetero_nonAS -o hetero_nonAS_stat.txt
     SNVAS stat -i PEsample.vcf -b PEsample_peaks_sorted.bed -t homo -o homo_stat.txt


Single-end ChIP-seq dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Simply run::

     cd SE_demo
     SNVAS call -b SEsample_peaks_sorted.bed -t SEsample_peaks_sorted.bam -c SEcontrol_peaks_sorted.bam -o SEsample.vcf

2. Filtering results::

     SNVAS filter -i SEsample.vcf -o hetero.vcf -t hetero -d 20 -q 100
     SNVAS filter -i SEsample.vcf -o hetero_AS.vcf -t hetero_AS -d 20 -q 100
     SNVAS filter -i SEsample.vcf -o hetero_nonAS.vcf -t hetero_nonAS -d 20 -q 100
     SNVAS filter -i SEsample.vcf -o homo.vcf -t homo -d 20 -q 10

3. Summarize predicted snvs number per kb and ts/tv ratio for different GQ cutoff::

     SNVAS stat -i SEsample.vcf -b SEsample_peaks_sorted.bed -t hetero -o hetero_stat.txt
     SNVAS stat -i SEsample.vcf -b SEsample_peaks_sorted.bed -t hetero_AS -o hetero_AS_stat.txt
     SNVAS stat -i SEsample.vcf -b SEsample_peaks_sorted.bed -t hetero_nonAS -o hetero_nonAS_stat.txt
     SNVAS stat -i SEsample.vcf -b SEsample_peaks_sorted.bed -t homo -o homo_stat.txt
