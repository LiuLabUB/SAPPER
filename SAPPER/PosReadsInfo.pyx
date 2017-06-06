# Time-stamp: <2017-06-05 16:42:21 Tao Liu>

"""Module for SAPPER BAMParser class

Copyright (c) 2017 Tao Liu <tliu4@buffalo.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included
with the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: tliu4@buffalo.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import logging
import struct
from struct import unpack
import gzip
import io
from collections import Counter


from SAPPER.Constants import *
from SAPPER.Stat import CalModel_Homo, CalModel_Heter_noAS, CalModel_Heter_AS, calculate_GQ, calculate_GQ_heterASsig

from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t, float32_t


cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    void free(void *ptr)
    int strcmp(char *a, char *b)
    char * strcpy(char *a, char *b)
    long atol(char *bytes)
    int atoi(char *bytes)

# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <tliu4@buffalo.edu>"
__doc__ = "All Parser classes"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class PosReadsInfo:
    cdef:
        long ref_pos
        bytes ref_nt
        bytes alt_nt
        bool filterout          # if true, do not output

        dict bq_set_T     #{A:[], C:[], G:[], T:[], N:[]} for treatment
        dict bq_set_C
        dict n_reads_T    #{A:[], C:[], G:[], T:[], N:[]} for treatment
        dict n_reads_C
        dict n_reads

        bytes top1nt
        bytes top2nt
        float top12NT_ratio

        double lnL_homo_major,lnL_heter_AS,lnL_heter_noAS,lnL_homo_minor
        double BIC_homo_major,BIC_heter_AS,BIC_heter_noAS,BIC_homo_minor
        int heter_noAS_kc, heter_noAS_ki
        int heter_AS_kc, heter_AS_ki
        double heter_AS_alleleratio

        int GQ_homo_major,GQ_heter_noAS,GQ_heter_AS  #phred scale of prob by standard formular
        int GQ_heter_ASsig #phred scale of prob, to measure the difference between AS and noAS

        int GQ
        
        str GT
        str type

        bool hasfermiinfor #if no fermi bam overlap in the position, false; if fermi bam in the position GT: N, false; if anyone of top2NT is not in fermi GT NTs, false;
        bytearray fermiNTs # 

    def __cinit__ ( self ):
        self.bq_set_T = { b'A':[], b'C':[], b'G':[], b'T':[], b'N':[] }
        self.bq_set_C = { b'A':[], b'C':[], b'G':[], b'T':[], b'N':[] }
        self.n_reads_T = { b'A':0, b'C':0, b'G':0, b'T':0, b'N':0 }
        self.n_reads_C = { b'A':0, b'C':0, b'G':0, b'T':0, b'N':0 }
        self.n_reads = { b'A':0, b'C':0, b'G':0, b'T':0, b'N':0 }
        self.filterout = False
        self.GQ = 0
        self.GT = "unsure"
        self.alt_nt = b'.'
        
    def __init__ ( self, ref_pos, ref_nt ):
        self.ref_pos = ref_pos
        self.ref_nt = ref_nt

    cpdef filterflag ( self ):
        return self.filterout

    cpdef apply_GQ_cutoff ( self, int min_homo_GQ = 50, int min_heter_GQ = 100 ):
        if self.type.startswith('homo') and self.GQ < min_homo_GQ:
            self.filterout = True
        elif self.type.startswith('heter') and self.GQ < min_heter_GQ:
            self.filterout = True
        return

    cpdef add_T ( self, int read_index, bytes read_nt, int read_bq ):
        self.bq_set_T[ read_nt ].append( (read_index, read_bq) )
        self.n_reads_T[ read_nt ] += 1
        self.n_reads[ read_nt ] += 1

    cpdef add_C ( self, int read_index, bytes read_nt, int read_bq ):
        self.bq_set_C[ read_nt ].append( (read_index, read_bq) )
        self.n_reads_C[ read_nt ] += 1
        self.n_reads[ read_nt ] += 1

    cpdef raw_read_depth ( self ):
        return sum( self.n_reads.values() )

    cpdef update_top_nts ( self, float min_top12NT_ratio = 0.8 ):
        """Identify top1 and top2 NT.  the ratio of (top1+top2)/total
        """
        cdef:
            float r
        [self.top1nt, self.top2nt] = sorted(self.n_reads, key=self.n_reads.get, reverse=True)[:2]
        
        self.top12NT_ratio = self.n_reads[ self.top1nt ] + self.n_reads[ self.top2nt ] /  sum( self.n_reads.values() )
        if self.top12NT_ratio < min_top12NT_ratio:
            self.filterout = True
        return

    cpdef top12NT ( self ):
        print ( self.ref_pos, self.ref_nt)
        print ("Top1NT",self.top1nt, "Treatment", self.bq_set_T[self.top1nt], "Control", self.bq_set_C[self.top1nt])
        print ("Top2NT",self.top2nt, "Treatment", self.bq_set_T[self.top2nt], "Control", self.bq_set_C[self.top2nt])
    
    cpdef call_GT ( self ):
        """Require update_top_nts being called.
        """
        cdef:
            np.ndarray[np.int32_t, ndim=1] top1_bq_T
            np.ndarray[np.int32_t, ndim=1] top2_bq_T
            np.ndarray[np.int32_t, ndim=1] top1_bq_C
            np.ndarray[np.int32_t, ndim=1] top2_bq_C
            int i
            list top1_bq_T_l
            list top2_bq_T_l
            list top1_bq_C_l
            list top2_bq_C_l

        if self.filterout:
            return

        top1_bq_T_l = self.bq_set_T[ self.top1nt ] 
        top2_bq_T_l = self.bq_set_T[ self.top2nt ] 
        top1_bq_C_l = self.bq_set_C[ self.top1nt ] 
        top2_bq_C_l = self.bq_set_C[ self.top2nt ] 

        top1_bq_T = np.zeros( len(top1_bq_T_l),dtype="int32")
        top2_bq_T = np.zeros( len(top2_bq_T_l),dtype="int32")
        top1_bq_C = np.zeros( len(top1_bq_C_l),dtype="int32")
        top2_bq_C = np.zeros( len(top2_bq_C_l),dtype="int32")

        for i in range( len(top1_bq_T_l) ):
            top1_bq_T[ i ] = top1_bq_T_l[ i ][ 1 ]
        for i in range( len(top2_bq_T_l) ):
            top2_bq_T[ i ] = top2_bq_T_l[ i ][ 1 ]
        for i in range( len(top1_bq_C_l) ):
            top1_bq_C[ i ] = top1_bq_C_l[ i ][ 1 ]
        for i in range( len(top2_bq_C_l) ):
            top2_bq_C[ i ] = top2_bq_C_l[ i ][ 1 ]
            
        (self.lnL_homo_major, self.BIC_homo_major) = CalModel_Homo( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_homo_minor, self.BIC_homo_minor) = CalModel_Homo( top2_bq_T, top2_bq_C, top1_bq_T, top1_bq_C )
        (self.lnL_heter_noAS, self.BIC_heter_noAS) = CalModel_Heter_noAS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_heter_AS, self.BIC_heter_AS)     = CalModel_Heter_AS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )

        self.GQ_homo_major = 0
        self.GQ_heter_noAS = 0
        self.GQ_heter_AS = 0
        self.GQ_heter_ASsig = 0
        
        # assign GQ, GT, and type
        if self.ref_nt != self.top1nt and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
            self.type = "homo"
            self.GQ_homo_major = calculate_GQ( self.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS )
            self.GQ = self.GQ_homo_major
            self.GT = "1/1"
            self.alt_nt = self.top1nt
        elif self.BIC_heter_noAS < self.BIC_homo_major and self.BIC_heter_noAS < self.BIC_homo_minor and self.BIC_heter_noAS < self.BIC_heter_AS+1e-8 :
            self.type = "heter_noAS"
            self.GQ_heter_noAS= calculate_GQ( self.lnL_heter_noAS, self.lnL_homo_major, self.lnL_homo_minor)
            self.GQ = self.GQ_heter_noAS
        elif self.BIC_heter_AS < self.BIC_homo_major and self.BIC_heter_AS < self.BIC_homo_minor and self.BIC_heter_AS < self.BIC_heter_noAS:
            self.type = "heter_AS"
            self.GQ_heter_AS = calculate_GQ( self.lnL_heter_AS, self.lnL_homo_major, self.lnL_homo_minor)
            self.GQ_heter_ASsig = calculate_GQ_heterASsig( self.lnL_heter_AS, self.lnL_heter_noAS )
            self.GQ = self.GQ_heter_AS
        elif self.ref_nt == self.top1nt and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
            self.type = "homo_ref"
            # we do not calculate GQ if type is homo_ref
            self.GT = "0/0"
            self.filterout = True
        else:
            self.type="unsure"
            self.filterout = True

        if self.type.startswith( "heter" ):
            if self.ref_nt == self.top1nt:
                self.alt_nt = self.top2nt
                self.GT = "0/1"
            elif self.ref_nt == self.top2nt:
                self.alt_nt = self.top1nt
                self.GT = "0/1"
            else:
                self.alt_nt = self.top1nt+b','+self.top2nt
                self.GT = "1/2"

        return

    cpdef to_vcf ( self ):
        """Output REF,ALT,QUAL,FILTER,INFO,FORMAT, SAMPLE columns.
        """
        cdef:
            str vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample

        # T       C       .       .       MinBIC_model=heter_AS;raw_depth_ChIP=56;raw_depth_input=0;DP_ChIP=46;DP_input=0;fermiNTs=CT;top1=33T;top2=13C;top1input=0T;top2input=0C;top1raw=40T;top2raw=16C;top1inputraw=0T;top2inputraw=0C;lnL_homo_major=-117.897;lnL_homo_minor=-296.115;lnL_heter_noAS=-31.7978;lnL_heter_AS=-29.4305;BIC_homo_major=235.794;BIC_homo_minor=592.229;BIC_heter_noAS=67.4243;BIC_heter_AS=66.5183;GQ_homo=0;GQ_heter_noAS=0;GQ_heter_AS=384;GQ_heter_ASsig=10;Allele_ratio_heter_AS=0.717391      GT:DP:GQ        0|1:46:384

        vcf_ref = self.ref_nt.decode()
        vcf_alt = self.alt_nt.decode()
        vcf_qual = "%d" % self.GQ
        vcf_filter = "."
        vcf_info = (b"MinBIC_model=%s;DP_ChIP=%d;DP_input=%d;fermiNTs=%d;top1=%d%s;top2=%d%s;top1input=%d%s;top2input=%d%s;lnL_homo_major=%f;lnL_homo_minor=%f;lnL_heter_noAS=%f;lnL_heter_AS=%f;BIC_homo_major=%f;BIC_homo_minor=%f;BIC_heter_noAS=%f;BIC_heter_AS=%f;GQ_homo=%d;GQ_heter_noAS=%d;GQ_heter_AS=%d;GQ_heter_ASsig=%d;Allele_ratio_heter_AS=%f" % \
            (self.type.encode(), sum( self.n_reads_T.values() ), sum( self.n_reads_C.values() ), 0, 
             self.n_reads_T[self.top1nt], self.top1nt, self.n_reads_T[self.top2nt], self.top2nt,
             self.n_reads_C[self.top1nt], self.top1nt, self.n_reads_C[self.top2nt], self.top2nt,
             self.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS, self.lnL_heter_AS,
             self.BIC_homo_major, self.BIC_homo_minor, self.BIC_heter_noAS, self.BIC_heter_AS,
             self.GQ_homo_major, self.GQ_heter_noAS, self.GQ_heter_AS, self.GQ_heter_ASsig, self.n_reads_T[self.top1nt]/(self.n_reads_T[self.top1nt]+self.n_reads_T[self.top2nt])
             )).decode()
        vcf_format = "GT:DP:GQ"
        vcf_sample = "%s:%d:%d" % (self.GT, self.raw_read_depth(), self.GQ )
        return "\t".join( ( vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample ) )
