# Time-stamp: <2017-06-09 16:52:06 Tao Liu>

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
import sys
from functools import partial
import multiprocessing as mp

from time import time
from math import ceil

# ------------------------------------
# own python modules
# ------------------------------------
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
    peakbedfile = args.peakbed
    peaktfile = args.tfile[0]
    peakcfile = args.cfile[0]
    top2ntminr = args.top2ntMinRatio
    NP = args.np

    # parameter for assembly
    fermiOverlapMinRatio = args.fermiOverlapMinRatio
    fermiOff = args.fermiOff
    
    peakio = open( peakbedfile )
    peaks = PeakIO()
    i = 0
    for l in peakio:
        fs = l.rstrip().split()
        i += 1
        peaks.add( fs[0].encode(), int(fs[1]), int(fs[2]), name=b"%d" % i )
    peaks.sort()

    chrs = peaks.get_chr_names()

    tbam = BAMParser( peaktfile )
    cbam = BAMParser( peakcfile )

    assert tbam.get_chromosomes() == cbam.get_chromosomes(), Exception("Treatment and Control BAM files have different orders of sorted chromosomes! Please check BAM Headers and re-sort BAM files.")

    ra_collections = []

    ovcf = open(args.ofile, "w")
    
    ovcf.write ( VCFHEADER % (datetime.date.today().strftime("%Y%m%d"), SAPPER_VERSION, " ".join(sys.argv[1:]) ) + "\n" )
    ovcf.write ( "\t".join( ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE") ) + "\n" )

    # to get time
    t_get_ra = 0
    t_get_pri = 0
    t_call_top2nt = 0
    t_call_lnL = 0
    t_call_GQ = 0
    t_call_to_vcf = 0
    t0 = time()

    for chrom in tbam.get_chromosomes():
        peaks_chr = peaks.get_data_from_chrom( chrom )
        for peak in peaks_chr:
            # note, when we extract reads from BAM within a peak
            # region, we assume BAM should be sorted and the BAM
            # should be generated from "samtools view -L" process.
            #t0 = time()
            ra_collection = RACollection( chrom, peak, tbam.get_reads_in_region( chrom, peak["start"], peak["end"] ), cbam.get_reads_in_region( chrom, peak["start"], peak["end"]) )
            ra_collection.remove_outliers( percent = 5 )
            if not fermiOff:
                # invoke fermi to assemble local sequence and filter out those can not be mapped to unitigs.
                unitigs = ra_collection.fermi_assemble( fermiOverlapMinRatio )
                ra_collection.filter_RAs_w_unitigs( unitigs )
            
            s = ra_collection.get_peak_REFSEQ()
            #t_get_ra += time() - t0
            # multiprocessing the following part

            if NP > 1:
                # divide right-left into NP parts
                window_size = ceil( ( ra_collection["right"] - ra_collection["left"] ) / NP )

                P = mp.Pool( NP )

                # this partial function will only be used in multiprocessing
                p_call_variants_at_range =  partial(call_variants_at_range, chrom=chrom, s=s, ra_collection=ra_collection, top2ntminr=top2ntminr)

                ranges = []
                
                for i in range( NP ):
                    l = i * window_size + ra_collection["left"]
                    r = min( (i + 1) * window_size + ra_collection["left"], ra_collection["right"] )
                    ranges.append( (l, r) )

                results = P.map( p_call_variants_at_range, ranges )
                P.close()
                P.join()
                for i in range( NP ):
                    ovcf.write( results[ i ] )
            else:
                for i in range( ra_collection["left"], ra_collection["right"] ):
                    ref_nt = chr(s[ i-ra_collection["left"] ] ).encode()

                    #t0 = time()
                    PRI = ra_collection.get_PosReadsInfo_ref_pos ( i, ref_nt )
                    #t_get_pri += time() - t0

                    #if PRI.raw_read_depth() == 0: # skip if coverage is 0
                    #    return None
                    #t0 = time()
                    #PRI.update_top_nts( top2ntminr )
                    #t_call_top2nt += time() - t0

                    #t0 = time()
                    #PRI.compute_lnL()
                    #t_call_lnL += time() - t0

                    #t0 = time()
                    #PRI.compute_GQ()
                    #t_call_GQ += time() - t0

                    #if not PRI.filterflag():
                    #    #t0 = time()
                    #    result=PRI.to_vcf()
                    #    #t_call_to_vcf += time() - t0
                    result  = call_variants_at_given_pos( PRI, top2ntminr )
                    #t_call += time() - t0
                    if result:
                        ovcf.write( "\t".join( ( chrom.decode(), str(i+1), ".", result ) ) + "\n" )

    #print ("time to read RAcollection from BAM:",t_get_ra)
    #print ("time to get reads information for each position:",t_get_pri)
    #print ("time to call top2 NTs:",t_call_top2nt)
    #print ("time to compute lnL:",t_call_lnL)
    #print ("time to compute GQ:",t_call_GQ)
    #print ("time to convert to vcf:",t_call_to_vcf)
    return

def call_variants_at_range ( lr, chrom, s, ra_collection, top2ntminr ):
    result = ""
    for i in range( lr[ 0 ], lr[ 1 ] ):
        ref_nt = chr(s[ i-ra_collection["left"] ] ).encode()
        PRI = ra_collection.get_PosReadsInfo_ref_pos ( i, ref_nt ) 
        if PRI.raw_read_depth() == 0: # skip if coverage is 0
            continue
        PRI.update_top_nts( top2ntminr )
        PRI.call_GT()
        PRI.apply_GQ_cutoff()
        if not PRI.filterflag():
            result += "\t".join( ( chrom.decode(), str(i+1), ".", PRI.to_vcf() ) ) + "\n"
    return result

def call_variants_at_given_pos ( PRI, top2ntminr ):
    if PRI.raw_read_depth() == 0: # skip if coverage is 0
        return None
    PRI.update_top_nts( top2ntminr )
    PRI.call_GT()
    PRI.apply_GQ_cutoff()
    if not PRI.filterflag():
        return PRI.to_vcf()
    return None
