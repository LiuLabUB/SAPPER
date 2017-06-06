# Time-stamp: <2017-06-05 16:31:17 Tao Liu>

"""Description: sapper call

Copyright (c) 2017 Tao Liu <tliu4@buffalo.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status: release candidate
@version: $Id$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import logging
import datetime
import tempfile
import sys

# ------------------------------------
# own python modules
# ------------------------------------
#from MACS2.OptValidator import opt_validate
#from MACS2.OutputWriter import *
#from MACS2.Prob import binomial_cdf_inv
#from MACS2.PeakModel import PeakModel,NotEnoughPairsException
#from MACS2.PeakDetect import PeakDetect
#from MACS2.Constants import *
from SAPPER.Constants import *
from SAPPER.PeakIO import PeakIO
from SAPPER.BAM import BAMParser
from SAPPER.RACollection import RACollection


VCFHEADER="""##fileformat=VCFv4.1
##fileDate=%s
##source=SAPPER_V%s
##Program_Args=%s
##INFO=<ID=MinBIC_model,Number=.,Type=String,Description="Model with minimum BIC value">
##INFO=<ID=DP_ChIP,Number=1,Type=Integer,Description="Read depth in ChIP-seq data">
##INFO=<ID=DP_input,Number=1,Type=Integer,Description="Read depth in input data">
##INFO=<ID=fermiNTs,Number=.,Type=String,Description="Nucleotides from the genotype information of fermi assembly result">
##INFO=<ID=top1,Number=.,Type=String,Description="Read depth of top1 nucleotide in ChIP-seq data">
##INFO=<ID=top2,Number=.,Type=String,Description="Read depth of top2 nucleotide in ChIP-seq data">
##INFO=<ID=top1input,Number=.,Type=String,Description="Read depth of top1 nucleotide in input data">
##INFO=<ID=top2input,Number=.,Type=String,Description="Read depth of top2 nucleotide in input data">
##INFO=<ID=lnL_homo_major,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with major allele model">
##INFO=<ID=lnL_homo_minor,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with minor allele model">
##INFO=<ID=lnL_heter_noAS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with no allele-specific model">
##INFO=<ID=lnL_heter_AS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with allele-specific model">
##INFO=<ID=BIC_homo_major,Number=1,Type=Float,Description="BIC value of homozygous with major allele model">
##INFO=<ID=BIC_homo_minor,Number=1,Type=Float,Description="BIC value of homozygous with minor allele model">
##INFO=<ID=BIC_heter_noAS,Number=1,Type=Float,Description="BIC value of heterozygous with no allele-specific model">
##INFO=<ID=BIC_heter_AS,Number=1,Type=Float,Description="BIC value of heterozygous with allele-specific model">
##INFO=<ID=GQ_homo,Number=1,Type=Integer,Description="Genotype quality of homozygous with major allele model">
##INFO=<ID=GQ_heter_noAS,Number=1,Type=Integer,Description="Genotype quality of heterozygous with no allele-specific model">
##INFO=<ID=GQ_heter_AS,Number=1,Type=Integer,Description="Genotype quality of heterozygous with allele-specific model">
##INFO=<ID=GQ_heter_ASsig,Number=1,Type=Integer,Description="Genotype quality of allele-specific significance compared with no allele-specific model">
##INFO=<ID=Allele_ratio_heter_AS,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="read depth after filtering">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality score">"""

# ------------------------------------
# Main function
# ------------------------------------

def check_names(treat, control, error_stream):
    """check common chromosome names"""
    tchrnames = set(treat.get_chr_names())
    cchrnames = set(control.get_chr_names())
    commonnames = tchrnames.intersection(cchrnames)
    if len(commonnames)==0:
        error_stream("No common chromosome names can be found from treatment and control! Check your input files! MACS will quit...")
        error_stream("Chromosome names in treatment: %s" % ",".join(sorted(tchrnames)))
        error_stream("Chromosome names in control: %s" % ",".join(sorted(cchrnames)))
        sys.exit()


def run( args ):
    """The Main function/pipeline for SAPPER

    """
    # Parse options...
    #options = opt_validate( args )
    # end of parsing commandline options
    
    peakbedfile = args.peakbed
    peaktfile = args.tfile[0]
    peakcfile = args.cfile[0]

    peakio = open( peakbedfile )
    peaks = PeakIO()
    i = 0
    for l in peakio:
        fs = l.rstrip().split()
        i += 1
        peaks.add( fs[0].encode(), int(fs[1]), int(fs[2]), name=b"%d" % i )
    peaks.sort()

    chrs = peaks.get_chr_names()

    b = BAMParser( peaktfile )
    c = BAMParser( peakcfile )

    ra_collections = []

    o = open(args.ofile, "w")
    
    o.write ( VCFHEADER % (datetime.date.today().strftime("%Y%m%d"), SAPPER_VERSION, " ".join(sys.argv[1:]) ) + "\n" )
    o.write ( "\t".join( ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE") ) + "\n" )

    for chrom in chrs:
        peaks_chr = peaks.get_data_from_chrom( chrom )
        for peak in peaks_chr:
            # note, when we extract reads from BAM within a peak
            # region, we assume BAM should be sorted and the BAM
            # should be generated from "samtools view -L" process.
            ra_collection = RACollection( chrom, peak, b.get_reads_in_region( chrom, peak["start"], peak["end"] ), c.get_reads_in_region( chrom, peak["start"], peak["end"]) )
            ra_collection.remove_outliers( percent = 5 )
            s = ra_collection.get_peak_REFSEQ()
            for i in range( ra_collection["left"], ra_collection["right"] ):
                ref_nt = chr(s[ i-ra_collection["left"] ] ).encode()
                PRI = ra_collection.get_PosReadsInfo_ref_pos ( i, ref_nt ) 
                PRI.update_top_nts()
                PRI.call_GT()
                PRI.apply_GQ_cutoff()
                if not PRI.filterflag():
                    o.write ( "\t".join( ( chrom.decode(), str(i+1), ".", PRI.to_vcf() ) ) + "\n" )
    return
