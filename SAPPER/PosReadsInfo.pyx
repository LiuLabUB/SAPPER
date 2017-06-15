# Time-stamp: <2017-06-15 17:24:06 Tao Liu>

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
        bytes ref_allele
        bytes alt_allele
        bool filterout          # if true, do not output

        dict bq_set_T     #{A:[], C:[], G:[], T:[], N:[]} for treatment
        dict bq_set_C
        dict n_reads_T    #{A:[], C:[], G:[], T:[], N:[]} for treatment
        dict n_reads_C
        dict n_reads

        bytes top1allele
        bytes top2allele
        float top12alleles_ratio

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
        str mutation_type       # SNV or Insertion or Deletion

        bool hasfermiinfor #if no fermi bam overlap in the position, false; if fermi bam in the position GT: N, false; if anyone of top2allele is not in fermi GT NTs, false;
        bytearray fermiNTs # 

    def __cinit__ ( self ):
        self.filterout = False
        self.GQ = 0
        self.GT = "unsure"
        self.alt_allele = b'.'
        
    def __init__ ( self, long ref_pos, bytes ref_allele ):
        self.ref_pos = ref_pos
        self.ref_allele = ref_allele
        self.bq_set_T = { ref_allele:[],b'A':[], b'C':[], b'G':[], b'T':[], b'N':[] }
        self.bq_set_C = { ref_allele:[],b'A':[], b'C':[], b'G':[], b'T':[], b'N':[] }
        self.n_reads_T = { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0 }
        self.n_reads_C = { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0 }
        self.n_reads =  { ref_allele:0,b'A':0, b'C':0, b'G':0, b'T':0, b'N':0 }

    def __getstate__ ( self ):
        return ( self.ref_pos, self.ref_allele, self.alt_allele, self.filterout,
                 self.bq_set_T, self.bq_set_C, self.n_reads_T, self.n_reads_C, self.n_reads,
                 self.top1allele, self.top2allele, self.top12alleles_ratio,
                 self.lnL_homo_major, self.lnL_heter_AS, self.lnL_heter_noAS, self.lnL_homo_minor,
                 self.BIC_homo_major, self.BIC_heter_AS, self.BIC_heter_noAS, self.BIC_homo_minor,
                 self.heter_noAS_kc, self.heter_noAS_ki,
                 self.heter_AS_kc, self.heter_AS_ki,
                 self.heter_AS_alleleratio,
                 self.GQ_homo_major, self.GQ_heter_noAS, self.GQ_heter_AS,
                 self.GQ_heter_ASsig,
                 self.GQ,
                 self.GT,
                 self.type,
                 self.hasfermiinfor,
                 self.fermiNTs )

    def __setstate__ ( self, state ):
        ( self.ref_pos, self.ref_allele, self.alt_allele, self.filterout,
          self.bq_set_T, self.bq_set_C, self.n_reads_T, self.n_reads_C, self.n_reads,
          self.top1allele, self.top2allele, self.top12alleles_ratio,
          self.lnL_homo_major, self.lnL_heter_AS, self.lnL_heter_noAS, self.lnL_homo_minor,
          self.BIC_homo_major, self.BIC_heter_AS, self.BIC_heter_noAS, self.BIC_homo_minor,
          self.heter_noAS_kc, self.heter_noAS_ki,
          self.heter_AS_kc, self.heter_AS_ki,
          self.heter_AS_alleleratio,
          self.GQ_homo_major, self.GQ_heter_noAS, self.GQ_heter_AS,
          self.GQ_heter_ASsig,
          self.GQ,
          self.GT,
          self.type,
          self.hasfermiinfor,
          self.fermiNTs ) = state

    cpdef filterflag ( self ):
        return self.filterout

    cpdef apply_GQ_cutoff ( self, int min_homo_GQ = 50, int min_heter_GQ = 100 ):
        if self.filterout:
            return
        if self.type.startswith('homo') and self.GQ < min_homo_GQ:
            self.filterout = True
        elif self.type.startswith('heter') and self.GQ < min_heter_GQ:
            self.filterout = True
        return

    cpdef add_T ( self, int read_index, bytes read_allele, int read_bq ):
        if not self.bq_set_T.has_key( read_allele ):
            self.bq_set_T[read_allele] = []
            self.bq_set_C[read_allele] = []
            self.n_reads_T[read_allele] = 0
            self.n_reads_C[read_allele] = 0
            self.n_reads[read_allele] = 0
        self.bq_set_T[read_allele].append( read_bq )
        self.n_reads_T[ read_allele ] += 1
        self.n_reads[ read_allele ] += 1

    cpdef add_C ( self, int read_index, bytes read_allele, int read_bq ):
        if not self.bq_set_C.has_key( read_allele ):
            self.bq_set_T[read_allele] = []
            self.bq_set_C[read_allele] = []
            self.n_reads_T[read_allele] = 0
            self.n_reads_C[read_allele] = 0
            self.n_reads[read_allele] = 0
        self.bq_set_C[read_allele].append( read_bq )
        self.n_reads_C[ read_allele ] += 1
        self.n_reads[ read_allele ] += 1

    cpdef raw_read_depth ( self ):
        return sum( self.n_reads.values() )

    cpdef update_top_alleles ( self, float min_top12alleles_ratio = 0.8 ):
        """Identify top1 and top2 NT.  the ratio of (top1+top2)/total
        """
        cdef:
            float r
        [self.top1allele, self.top2allele] = sorted(self.n_reads, key=self.n_reads.get, reverse=True)[:2]
        
        self.top12alleles_ratio = ( self.n_reads[ self.top1allele ] + self.n_reads[ self.top2allele ] ) /  sum( self.n_reads.values() )
        if self.top12alleles_ratio < min_top12alleles_ratio:
            self.filterout = True
        if self.top1allele == self.ref_allele and self.n_reads[ self.top2allele ] == 0:
            # This means this position only contains top1allele which is the ref_allele. So the GT must be 0/0
            self.filterout = True
        return

    cpdef top12alleles ( self ):
        print ( self.ref_pos, self.ref_allele)
        print ("Top1allele",self.top1allele, "Treatment", self.bq_set_T[self.top1allele], "Control", self.bq_set_C[self.top1allele])
        print ("Top2allele",self.top2allele, "Treatment", self.bq_set_T[self.top2allele], "Control", self.bq_set_C[self.top2allele])
    
    cpdef call_GT ( self ):
        """Require update_top_alleles being called.
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
            list tmp_mutation_type
            bytes tmp_alt

        if self.filterout:
            return

        #top1_bq_T_l = self.bq_set_T[ self.top1allele ] 
        #top2_bq_T_l = self.bq_set_T[ self.top2allele ] 
        #top1_bq_C_l = self.bq_set_C[ self.top1allele ] 
        #top2_bq_C_l = self.bq_set_C[ self.top2allele ] 

        #top1_bq_T = np.zeros( len(top1_bq_T_l),dtype="int32")
        #top2_bq_T = np.zeros( len(top2_bq_T_l),dtype="int32")
        #top1_bq_C = np.zeros( len(top1_bq_C_l),dtype="int32")
        #top2_bq_C = np.zeros( len(top2_bq_C_l),dtype="int32")
        #print self.top1allele
        #print self.top2allele
        top1_bq_T = np.array( self.bq_set_T[ self.top1allele ], dtype="int32" )
        top2_bq_T = np.array( self.bq_set_T[ self.top2allele ], dtype="int32" )
        top1_bq_C = np.array( self.bq_set_C[ self.top1allele ], dtype="int32" )
        top2_bq_C = np.array( self.bq_set_C[ self.top2allele ], dtype="int32" )

        #for i in range( len(top1_bq_T_l) ):
        #    top1_bq_T[ i ] = top1_bq_T_l[ i ][ 1 ]
        #for i in range( len(top2_bq_T_l) ):
        #    top2_bq_T[ i ] = top2_bq_T_l[ i ][ 1 ]
        #for i in range( len(top1_bq_C_l) ):
        #    top1_bq_C[ i ] = top1_bq_C_l[ i ][ 1 ]
        #for i in range( len(top2_bq_C_l) ):
        #    top2_bq_C[ i ] = top2_bq_C_l[ i ][ 1 ]
            
        (self.lnL_homo_major, self.BIC_homo_major) = CalModel_Homo( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_homo_minor, self.BIC_homo_minor) = CalModel_Homo( top2_bq_T, top2_bq_C, top1_bq_T, top1_bq_C )
        (self.lnL_heter_noAS, self.BIC_heter_noAS) = CalModel_Heter_noAS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_heter_AS, self.BIC_heter_AS)     = CalModel_Heter_AS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )

        self.GQ_homo_major = 0
        self.GQ_heter_noAS = 0
        self.GQ_heter_AS = 0
        self.GQ_heter_ASsig = 0
        
        # assign GQ, GT, and type
        if self.ref_allele != self.top1allele and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
            self.type = "homo"
            self.GQ_homo_major = calculate_GQ( self.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS )
            self.GQ = self.GQ_homo_major
            self.GT = "1/1"
            self.alt_allele = self.top1allele
        elif self.BIC_heter_noAS < self.BIC_homo_major and self.BIC_heter_noAS < self.BIC_homo_minor and self.BIC_heter_noAS < self.BIC_heter_AS+1e-8 :
            self.type = "heter_noAS"
            self.GQ_heter_noAS= calculate_GQ( self.lnL_heter_noAS, self.lnL_homo_major, self.lnL_homo_minor)
            self.GQ = self.GQ_heter_noAS
        elif self.BIC_heter_AS < self.BIC_homo_major and self.BIC_heter_AS < self.BIC_homo_minor and self.BIC_heter_AS < self.BIC_heter_noAS:
            self.type = "heter_AS"
            self.GQ_heter_AS = calculate_GQ( self.lnL_heter_AS, self.lnL_homo_major, self.lnL_homo_minor)
            self.GQ_heter_ASsig = calculate_GQ_heterASsig( self.lnL_heter_AS, self.lnL_heter_noAS )
            self.GQ = self.GQ_heter_AS
        elif self.ref_allele == self.top1allele and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
            self.type = "homo_ref"
            # we do not calculate GQ if type is homo_ref
            self.GT = "0/0"
            self.filterout = True
        else:
            self.type="unsure"
            self.filterout = True

        if self.type.startswith( "heter" ):
            if self.ref_allele == self.top1allele:
                self.alt_allele = self.top2allele
                self.GT = "0/1"
            elif self.ref_allele == self.top2allele:
                self.alt_allele = self.top1allele
                self.GT = "0/1"
            else:
                self.alt_allele = self.top1allele+b','+self.top2allele
                self.GT = "1/2"

        tmp_mutation_type = []
        for tmp_alt in self.alt_allele.split(b','):
            if tmp_alt == b'*':
                tmp_mutation_type.append( "Deletion" )
            elif len( tmp_alt ) > 1:
                tmp_mutation_type.append( "Insertion" )
            else:
                tmp_mutation_type.append( "SNV" )
        self.mutation_type = ",".join( tmp_mutation_type )


        return

    cpdef compute_lnL ( self ):
        """Require update_top_alleles being called.
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

        top1_bq_T = np.array( self.bq_set_T[ self.top1allele ], dtype="int32" )
        top2_bq_T = np.array( self.bq_set_T[ self.top2allele ], dtype="int32" )
        top1_bq_C = np.array( self.bq_set_C[ self.top1allele ], dtype="int32" )
        top2_bq_C = np.array( self.bq_set_C[ self.top2allele ], dtype="int32" )

        (self.lnL_homo_major, self.BIC_homo_major) = CalModel_Homo( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_homo_minor, self.BIC_homo_minor) = CalModel_Homo( top2_bq_T, top2_bq_C, top1_bq_T, top1_bq_C )
        (self.lnL_heter_noAS, self.BIC_heter_noAS) = CalModel_Heter_noAS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )
        (self.lnL_heter_AS, self.BIC_heter_AS)     = CalModel_Heter_AS( top1_bq_T, top1_bq_C, top2_bq_T, top2_bq_C )

        return

    cpdef compute_GQ ( self ):
        cdef:
            list tmp_mutation_type
            bytes tmp_alt

        if self.filterout:
            return
        self.GQ_homo_major = 0
        self.GQ_heter_noAS = 0
        self.GQ_heter_AS = 0
        self.GQ_heter_ASsig = 0
        
        # assign GQ, GT, and type
        if self.ref_allele != self.top1allele and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
            self.type = "homo"
            self.GQ_homo_major = calculate_GQ( self.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS )
            self.GQ = self.GQ_homo_major
            self.GT = "1/1"
            self.alt_allele = self.top1allele
        elif self.BIC_heter_noAS < self.BIC_homo_major and self.BIC_heter_noAS < self.BIC_homo_minor and self.BIC_heter_noAS < self.BIC_heter_AS+1e-8 :
            self.type = "heter_noAS"
            self.GQ_heter_noAS= calculate_GQ( self.lnL_heter_noAS, self.lnL_homo_major, self.lnL_homo_minor)
            self.GQ = self.GQ_heter_noAS
        elif self.BIC_heter_AS < self.BIC_homo_major and self.BIC_heter_AS < self.BIC_homo_minor and self.BIC_heter_AS < self.BIC_heter_noAS:
            self.type = "heter_AS"
            self.GQ_heter_AS = calculate_GQ( self.lnL_heter_AS, self.lnL_homo_major, self.lnL_homo_minor)
            self.GQ_heter_ASsig = calculate_GQ_heterASsig( self.lnL_heter_AS, self.lnL_heter_noAS )
            self.GQ = self.GQ_heter_AS
        elif self.ref_allele == self.top1allele and self.BIC_homo_major < self.BIC_homo_minor and self.BIC_homo_major < self.BIC_heter_noAS and self.BIC_homo_major < self.BIC_heter_AS:
            self.type = "homo_ref"
            # we do not calculate GQ if type is homo_ref
            self.GT = "0/0"
            self.filterout = True
        else:
            self.type="unsure"
            self.filterout = True

        if self.type.startswith( "heter" ):
            if self.ref_allele == self.top1allele:
                self.alt_allele = self.top2allele
                self.GT = "0/1"
            elif self.ref_allele == self.top2allele:
                self.alt_allele = self.top1allele
                self.GT = "0/1"
            else:
                self.alt_allele = self.top1allele+b','+self.top2allele
                self.GT = "1/2"

        tmp_mutation_type = []
        for tmp_alt in self.alt_allele.split(b','):
            if tmp_alt == 42: # means '*'
                tmp_mutation_type.append( "Deletion" )
            elif len( tmp_alt ) > 1:
                tmp_mutation_type.append( "Insertion" )
            else:
                tmp_mutation_type.append( "SNV" )
        self.mutation_type = ",".join( tmp_mutation_type )
        return

    cpdef to_vcf ( self ):
        """Output REF,ALT,QUAL,FILTER,INFO,FORMAT, SAMPLE columns.
        """
        cdef:
            str vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample

        vcf_ref = self.ref_allele.decode()
        vcf_alt = self.alt_allele.decode()
        vcf_qual = "%d" % self.GQ
        vcf_filter = "."
#        vcf_info = (b"M=%s;MT=%s;DPT=%d;DPC=%d;DP1T=%d%s;DP2T=%d%s;DP1C=%d%s;DP2C=%d%s;lnLHOMOMAJOR=%.4f;lnLHOMOMINOR=%.4f;lnLHETERNOAS=%.4f;lnLHETERAS=%.4f;BICHOMOMAJOR=%.4f;BICHOMOMINOR=%.4f;BICHETERNOAS=%.4f;BICHETERAS=%.4f;GQHOMO=%d;GQHETERNOAS=%d;GQHETERAS=%d;GQHETERASsig=%d;AR=%.4f" % \
        vcf_info = (b"M=%s;MT=%s;DPT=%d;DPC=%d;DP1T=%d%s;DP2T=%d%s;DP1C=%d%s;DP2C=%d%s;GQHOMO=%d;GQHETERNOAS=%d;GQHETERAS=%d;GQHETERASsig=%d;AR=%.4f" % \
            (self.type.encode(), self.mutation_type.encode(), sum( self.n_reads_T.values() ), sum( self.n_reads_C.values() ), 
             self.n_reads_T[self.top1allele], self.top1allele, self.n_reads_T[self.top2allele], self.top2allele,
             self.n_reads_C[self.top1allele], self.top1allele, self.n_reads_C[self.top2allele], self.top2allele,
#             self.lnL_homo_major, self.lnL_homo_minor, self.lnL_heter_noAS, self.lnL_heter_AS,
#             self.BIC_homo_major, self.BIC_homo_minor, self.BIC_heter_noAS, self.BIC_heter_AS,
             self.GQ_homo_major, self.GQ_heter_noAS, self.GQ_heter_AS, self.GQ_heter_ASsig, self.n_reads_T[self.top1allele]/(self.n_reads_T[self.top1allele]+self.n_reads_T[self.top2allele])
             )).decode()
        vcf_format = "GT:DP:GQ"
        vcf_sample = "%s:%d:%d" % (self.GT, self.raw_read_depth(), self.GQ )
        return "\t".join( ( vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, vcf_sample ) )
