# Time-stamp: <2017-07-25 16:05:57 Tao Liu>

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


VCFHEADER_0="""##fileformat=VCFv4.1
##fileDate=%s
##source=SAPPER_V%s
##Program_Args=%s
##INFO=<ID=M,Number=.,Type=String,Description="SAPPER Model with minimum BIC value">
##INFO=<ID=MT,Number=.,Type=String,Description="Mutation type: SNV/Insertion/Deletion">
##INFO=<ID=DPT,Number=1,Type=Integer,Description="Depth Treatment: Read depth in ChIP-seq data">
##INFO=<ID=DPC,Number=1,Type=Integer,Description="Depth Control: Read depth in control data">
##INFO=<ID=DP1T,Number=.,Type=String,Description="Read depth of top1 allele in ChIP-seq data">
##INFO=<ID=DP2T,Number=.,Type=String,Description="Read depth of top2 allele in ChIP-seq data">
##INFO=<ID=DP1C,Number=.,Type=String,Description="Read depth of top1 allele in control data">
##INFO=<ID=DP2C,Number=.,Type=String,Description="Read depth of top2 allele in control data">
##INFO=<ID=lnLHOMOMAJOR,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with major allele model">
##INFO=<ID=lnLHOMOMINOR,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of homozygous with minor allele model">
##INFO=<ID=lnLHETERNOAS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with no allele-specific model">
##INFO=<ID=lnLHETERAS,Number=1,Type=Float,Description="Log(e) scaled genotype likelihoods of heterozygous with allele-specific model">
##INFO=<ID=BICHOMOMAJOR,Number=1,Type=Float,Description="BIC value of homozygous with major allele model">
##INFO=<ID=BICHOMOMINOR,Number=1,Type=Float,Description="BIC value of homozygous with minor allele model">
##INFO=<ID=BICHETERNOAS,Number=1,Type=Float,Description="BIC value of heterozygous with no allele-specific model">
##INFO=<ID=BICHETERAS,Number=1,Type=Float,Description="BIC value of heterozygous with allele-specific model">
##INFO=<ID=GQHOMO,Number=1,Type=Integer,Description="Genotype quality of homozygous with major allele model">
##INFO=<ID=GQHETERNOAS,Number=1,Type=Integer,Description="Genotype quality of heterozygous with no allele-specific model">
##INFO=<ID=GQHETERAS,Number=1,Type=Integer,Description="Genotype quality of heterozygous with allele-specific model">
##INFO=<ID=GQHETERASsig,Number=1,Type=Integer,Description="Genotype quality of allele-specific significance compared with no allele-specific model">
##INFO=<ID=AR,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="read depth after filtering">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality score">"""

VCFHEADER="""##fileformat=VCFv4.1
##fileDate=%s
##source=SAPPER_V%s
##Program_Args=%s
##INFO=<ID=M,Number=.,Type=String,Description="SAPPER Model with minimum BIC value">
##INFO=<ID=MT,Number=.,Type=String,Description="Mutation type: SNV/Insertion/Deletion">
##INFO=<ID=DPT,Number=1,Type=Integer,Description="Depth Treatment: Read depth in ChIP-seq data">
##INFO=<ID=DPC,Number=1,Type=Integer,Description="Depth Control: Read depth in control data">
##INFO=<ID=DP1T,Number=.,Type=String,Description="Read depth of top1 allele in ChIP-seq data">
##INFO=<ID=DP2T,Number=.,Type=String,Description="Read depth of top2 allele in ChIP-seq data">
##INFO=<ID=DP1C,Number=.,Type=String,Description="Read depth of top1 allele in control data">
##INFO=<ID=DP2C,Number=.,Type=String,Description="Read depth of top2 allele in control data">
##INFO=<ID=DBIC,Number=.,Type=Float,Description="Difference of BIC of selected model and alternative model">
##INFO=<ID=BICHOMOMAJOR,Number=1,Type=Integer,Description="BIC of homozygous with major allele model">
##INFO=<ID=BICHOMOMINOR,Number=1,Type=Integer,Description="BIC of homozygous with minor allele model">
##INFO=<ID=BICHETERNOAS,Number=1,Type=Integer,Description="BIC of heterozygous with no allele-specific model">
##INFO=<ID=BICHETERAS,Number=1,Type=Integer,Description="BIC of heterozygous with allele-specific model">
##INFO=<ID=AR,Number=1,Type=Float,Description="Estimated allele ratio of heterozygous with allele-specific model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="read depth after filtering">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality score">
##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Phred-scaled genotype likelihoods for 00, 01, 11 genotype">"""

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
    top2allelesminr = args.top2allelesMinRatio
    NP = args.np
    min_homo_GQ = args.GQCutoffHomo
    min_heter_GQ = args.GQCutoffHetero

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
    
    ovcf.write ( VCFHEADER % (datetime.date.today().strftime("%Y%m%d"), SAPPER_VERSION, " ".join(sys.argv[1:] + ["-g",str(min_heter_GQ),"-G",str(min_homo_GQ)]) ) + "\n" )
    for (chrom, chrlength) in tbam.get_rlengths().items():
        ovcf.write( "##contig=<ID=%s,length=%d,assembly=NA>\n" % ( chrom.decode(), chrlength ) )


    ovcf.write ( "\t".join( ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE") ) + "\n" )

    # to get time
    t_get_ra = 0
    t_get_pri = 0
    t_call_top2alleles = 0
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

            # print ( "---begin of peak---")
            print ( "Peak:", chrom.decode(), peak["start"], peak["end"])

            # t0 = time()
            ra_collection = RACollection( chrom, peak, tbam.get_reads_in_region( chrom, peak["start"], peak["end"] ), cbam.get_reads_in_region( chrom, peak["start"], peak["end"]) )
            ra_collection.remove_outliers( percent = 5 )

            # print ( "Reads in Peak:")
            # print ( ra_collection.get_FASTQ().decode() )

            if not fermiOff:
                # invoke fermi to assemble local sequence and filter out those can not be mapped to unitigs.

                print ( " Assemble using fermi-lite")
                unitigs = ra_collection.fermi_assemble( fermiOverlapMinRatio )
                if unitigs == []:
                    print ( " failed to assemble unitigs")
                    continue              #pass this peak if there is no unitig from assembler

                print ( " Re-align reads with unitigs")
                (unitig_alns, reference_alns) = ra_collection.align_unitig_to_REFSEQ( unitigs )
                unitig_collection = ra_collection.remap_RAs_w_unitigs( unitigs, (unitig_alns, reference_alns) )                

                # print ( " Assembled unitigs collection:")
                # print ( " ", unitig_collection["chrom"].decode(), unitig_collection["left"], unitig_collection["right"], unitig_collection["URAs_left"], unitig_collection["URAs_right"], unitig_collection["count"] )
                # for (i, ua) in enumerate(unitig_collection["URAs_list"]):
                #     print ("  Unitig",i)
                #     print ( "  ", ua["chrom"].decode(), ua["lpos"], ua["rpos"], ua["count"] )
                #     print ( "  unitig seq:", ua["seq"].decode() )
                #     print ( "  unitig aln:", ua["unitig_aln"].decode() )
                #     print ( "  refere aln:", ua["reference_aln"].decode() )

            # print ( "---end of peak---" )
            
            s = ra_collection["peak_refseq"]
            #t_get_ra += time() - t0
            # multiprocessing the following part

            print ( " Call varants")
            if NP > 1:
                # divide right-left into NP parts
                window_size = ceil( ( ra_collection["right"] - ra_collection["left"] ) / NP )

                P = mp.Pool( NP )

                # this partial function will only be used in multiprocessing
                if not fermiOff:
                    p_call_variants_at_range =  partial(call_variants_at_range, chrom=chrom, s=s, collection=unitig_collection, top2allelesminr=top2allelesminr, min_homo_GQ = min_homo_GQ, min_heter_GQ = min_heter_GQ)
                else:
                    p_call_variants_at_range =  partial(call_variants_at_range, chrom=chrom, s=s, collection=ra_collection, top2allelesminr=top2allelesminr, min_homo_GQ = min_homo_GQ, min_heter_GQ = min_heter_GQ)

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

                    # skip if ref_nt is N
                    if ref_nt == b'N':
                        continue

                    #t0 = time()
                    if not fermiOff:
                        PRI = unitig_collection.get_PosReadsInfo_ref_pos ( i, ref_nt )
                    else:
                        PRI = ra_collection.get_PosReadsInfo_ref_pos ( i, ref_nt )
                    if PRI.raw_read_depth() == 0: # skip if coverage is 0
                        continue
                    PRI.update_top_alleles( top2allelesminr )
                    PRI.call_GT()
                    PRI.apply_GQ_cutoff(min_homo_GQ, min_heter_GQ)
                    # if i == 66176154 or i == 66176155:
                    #    PRI.top12alleles()
                    #    print ( PRI.to_vcf() )
                    #    print ( PRI.filterflag() )
                    if not PRI.filterflag():
                        ovcf.write( "\t".join( ( chrom.decode(), str(i+1), ".", PRI.to_vcf() ) ) + "\n" )

    #print ("time to read RAcollection from BAM:",t_get_ra)
    #print ("time to get reads information for each position:",t_get_pri)
    #print ("time to call top2 NTs:",t_call_top2alleles)
    #print ("time to compute lnL:",t_call_lnL)
    #print ("time to compute GQ:",t_call_GQ)
    #print ("time to convert to vcf:",t_call_to_vcf)
    return

def call_variants_at_range ( lr, chrom, s, collection, top2allelesminr, min_homo_GQ, min_heter_GQ ):
    result = ""
    for i in range( lr[ 0 ], lr[ 1 ] ):
        ref_nt = chr(s[ i-collection["left"] ] ).encode()
        if ref_nt == b'N':
            continue

        PRI = collection.get_PosReadsInfo_ref_pos ( i, ref_nt ) 
        if PRI.raw_read_depth() == 0: # skip if coverage is 0
            continue
        PRI.update_top_alleles( top2allelesminr )
        PRI.call_GT()
        PRI.apply_GQ_cutoff(min_homo_GQ, min_heter_GQ)
        if not PRI.filterflag():
            result += "\t".join( ( chrom.decode(), str(i+1), ".", PRI.to_vcf() ) ) + "\n"
    return result

# def call_variants_at_given_pos ( PRI, top2allelesminr ):
#     if PRI.raw_read_depth() == 0: # skip if coverage is 0
#         return None
#     PRI.update_top_alleles( top2allelesminr )
#     PRI.call_GT()
#     PRI.apply_GQ_cutoff()
#     if not PRI.filterflag():
#         return PRI.to_vcf()
#     return None
