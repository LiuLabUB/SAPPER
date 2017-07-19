# Time-stamp: <2017-07-19 16:02:18 Tao Liu>

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
from collections import Counter
from operator import itemgetter
from copy import copy

from SAPPER.Constants import *
from SAPPER.ReadAlignment import ReadAlignment
from SAPPER.PosReadsInfo import PosReadsInfo
from SAPPER.PeakIO import PeakIO
from SAPPER.UnitigRACollection import UnitigRAs, UnitigCollection
#from SAPPER.Alignment import SWalign

from cpython cimport bool
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t, int64_t

from libc.stdlib cimport malloc, free, realloc

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)
    #void *malloc(size_t size)
    void *calloc(size_t n, size_t size)
    #void free(void *ptr)
    int strcmp(char *a, char *b)
    char * strcpy(char *a, char *b)
    long atol(char *bytes)
    int atoi(char *bytes)


# --- fermi-lite functions ---
#define MAG_F_AGGRESSIVE 0x20 // pop variant bubbles (not default)
#define MAG_F_POPOPEN    0x40 // aggressive tip trimming (default)
#define MAG_F_NO_SIMPL   0x80 // skip bubble simplification (default)

cdef extern from "fml.h":
    ctypedef struct bseq1_t:
        int32_t l_seq
        char *seq
        char *qual # NULL-terminated strings; length expected to match $l_seq

    ctypedef struct magopt_t:
        int flag, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bdiff, max_bvtx, min_merge_len, trim_len, trim_depth
        float min_dratio1, max_bcov, max_bfrac

    ctypedef struct fml_opt_t:
        int n_threads        # number of threads; don't use multi-threading for small data sets
        int ec_k             # k-mer length for error correction; 0 for auto estimate
        int min_cnt, max_cnt # both occ threshold in ec and tip threshold in cleaning lie in [min_cnt,max_cnt]
        int min_asm_ovlp     # min overlap length during assembly
        int min_merge_len    # during assembly, don't explicitly merge an overlap if shorter than this value
        magopt_t mag_opt     # graph cleaning options

    ctypedef struct fml_ovlp_t:
        uint32_t len_, from_, id_, to_ 
        #unit32_t from  # $from and $to: 0 meaning overlapping 5'-end; 1 overlapping 3'-end
        #uint32_t id
        #uint32_t to    # $id: unitig number

    ctypedef struct fml_utg_t:
        int32_t len      # length of sequence
        int32_t nsr      # number of supporting reads
        char *seq        # unitig sequence
        char *cov        # cov[i]-33 gives per-base coverage at i
        int n_ovlp[2]    # number of 5'-end [0] and 3'-end [1] overlaps
        fml_ovlp_t *ovlp # overlaps, of size n_ovlp[0]+n_ovlp[1]

    void fml_opt_init(fml_opt_t *opt)
    fml_utg_t* fml_assemble(const fml_opt_t *opt, int n_seqs, bseq1_t *seqs, int *n_utg)
    void fml_utg_destroy(int n_utg, fml_utg_t *utg)
    void fml_utg_print(int n_utgs, const fml_utg_t *utg)

# cdef extern from "kstring.h":
#     ctypedef struct kstring_t:
#         uint32_t l, m
#         char *s

# cdef extern from "priv.h":
#     int ksa_bwt64(unsigned char *T, int64_t n, int k)

# cdef extern from "mag.h":
#     ctypedef struct magopt_t:
#         int flag, max_arc, n_iter, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bvtx, min_merge_len
#         float min_dratio0, min_dratio1
#         float max_bcov, max_bfrac
#     ctypedef struct ku128_t:
#         uint64_t x, y
#     ctypedef struct ku64_v:
#         size_t n, m
#         uint64_t *a
#     ctypedef struct ku128_v:
#         size_t n, m
#         ku128_t *a
#     ctypedef struct magv_t:
#         int len, nsr            #length; number supporting reads
#         uint32_t max_len        # allocated seq/cov size
#         uint64_t k[2]           # bi-interval
#         ku128_v nei[2]          # neighbors
#         char *seq
#         char *cov               # sequence and coverage
#         void *ptr               # additional information
#     ctypedef struct magv_v:
#         size_t n, m
#         magv_t *a
#     ctypedef struct mag_t:
#         magv_v v
#         float rdist             # read distance
#         int min_ovlp            # minimum overlap seen from the graph
#         void *h

#     magopt_t *mag_init_opt()
#     void mag_g_clean(mag_t *g, const magopt_t *opt)
#     void mag_g_print(const mag_t *g)
#     void mag_g_destroy(mag_t *g)
#     void mag_v_write(const magv_t *p, kstring_t *out);
    
# cdef extern from "fermi.h":
#    int fm6_api_correct(int kmer, int64_t l, char *_seq, char *_qual)
#    mag_t *fm6_api_unitig(int min_match, int64_t l, char *seq)


# --- end of fermi-lite functions ---

# --- smith-waterman alignment functions ---

cdef extern from "swalign.h":
    ctypedef struct seq_pair_t:
        char *a
        unsigned int alen
        char *b
        unsigned int blen
    ctypedef struct align_t:
        seq_pair_t *seqs
        int start_a
        int start_b
        int end_a
        int end_b
        int matches
        double score
    align_t *smith_waterman(seq_pair_t *problem)
    void destroy_seq_pair(seq_pair_t *pair)
    
# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <tliu4@buffalo.edu>"
__doc__ = "All Parser classes"

__DNACOMPLEMENT__ = b'\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQRSAUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff' # A trans table to convert A to T, C to G, G to C, and T to A.

__CIGARCODE__ = "MIDNSHP=X"

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

cdef class RACollection:
    """A collection of ReadAlignment objects and the corresponding
    PeakIO.

    """
    cdef:
        bytes chrom
        object peak             # A PeakIO object
        list RAlists           # contain ReadAlignment lists for treatment (0) and control (1)
        long left               # left position of peak
        long right              # right position of peak
        long length             # length of peak
        long RAs_left           # left position of all RAs in the collection
        long RAs_right          # right position of all RAs in the collection
        bool sorted             # if sorted by lpos
        bytes peak_refseq       # reference sequence in peak region b/w left and right
        bytes peak_refseq_ext   # reference sequence in peak region with extension on both sides b/w RAs_left and RAs_right

    def __init__ ( self, chrom, peak, RAlist_T, RAlist_C=[] ):
        """Create RACollection object by taking:

        1. peak: a PeakIO object indicating the peak region.

        2. RAlist: a python list of ReadAlignment objects containing
        all the reads overlapping the peak region. If no RAlist_C
        given, it will be [].

        """
        self.chrom = chrom
        self.peak = peak
        self.RAlists = [ RAlist_T, RAlist_C ]
        self.left = peak["start"]
        self.right = peak["end"]
        self.length =  self.right - self.left
        self.RAs_left = RAlist_T[0]["lpos"] # initial assignment of RAs_left
        self.RAs_right = RAlist_T[-1]["rpos"] # initial assignment of RAs_right
        self.sort()                           # it will set self.sorted = True
        # check RAs_left and RAs_right
        for ra in RAlist_T:
            if ra[ "lpos" ] < self.RAs_left:
                self.RAs_left = ra[ "lpos" ]
            if ra[ "rpos" ] > self.RAs_right:
                self.RAs_right = ra[ "rpos" ]

        for ra in RAlist_C:
            if ra[ "lpos" ] < self.RAs_left:
                self.RAs_left = ra[ "lpos" ]
            if ra[ "rpos" ] > self.RAs_right:
                self.RAs_right = ra[ "rpos" ]
        (self.peak_refseq, self.peak_refseq_ext) = self.__get_peak_REFSEQ()

    def __getitem__ ( self, keyname ):
        if keyname == "chrom":
            return self.chrom
        elif keyname == "left":
            return self.left
        elif keyname == "right":
            return self.right
        elif keyname == "RAs_left":
            return self.RAs_left
        elif keyname == "RAs_right":
            return self.RAs_right	    
        elif keyname == "length":
            return self.length
        elif keyname == "count":
            return len( self.RAlists[ 0 ] )+ len( self.RAlists[ 1 ] )
        elif keyname == "count_T":
            return len( self.RAlists[ 0 ] )
        elif keyname == "count_C":
            return len( self.RAlists[ 1 ] )
        elif keyname == "peak_refseq":
            return self.peak_refseq
        elif keyname == "peak_refseq_ext":
            return self.peak_refseq_ext
        else:
            raise KeyError("Unavailable key:", keyname)

    def __getstate__ ( self ):
        #return {"chrom":self.chrom, "peak":self.peak, "RAlists":self.RAlists,
        #        "left":self.left, "right":self.right, "length": self.length,
        #        "RAs_left":self.RAs_left, "RAs_right":self.RAs_right}
        return (self.chrom, self.peak, self.RAlists, self.left, self.right, self.length, self.RAs_left, self.RAs_right, self.peak_refseq, self.peak_refseq_ext)
        
    def __setstate__ ( self, state ):
        (self.chrom, self.peak, self.RAlists, self.left, self.right, self.length, self.RAs_left, self.RAs_right, self.peak_refseq, self.peak_refseq_ext) = state
        
    cpdef sort ( self ):
        """Sort RAs according to lpos. Should be used after realignment.

        """
        self.RAlists[ 0 ].sort(key=itemgetter("lpos"))
        self.RAlists[ 1 ].sort(key=itemgetter("lpos"))
        self.sorted = True
        return
        
    cpdef remove_outliers ( self, int percent = 5 ):
        """ Remove outliers with too many n_edits. The outliers with
        n_edits in top p% will be removed.

        Default: remove top 5% of reads that have too many differences
        with reference genome.
        """
        cdef:
            list n_edits_list
            object read         # ReadAlignment object
            int highest_n_edits
            list new_RAlist
            int i

        n_edits_list = []
        for ralist in self.RAlists:
            for read in ralist:
                n_edits_list.append( read["n_edits"] )
        n_edits_list.sort()
        highest_n_edits = n_edits_list[ int( len( n_edits_list ) * (1 - percent * .01) ) ]

        for i in ( range(len(self.RAlists)) ):
            new_RAlist = []
            for read in self.RAlists[ i ]:
                if read["n_edits"] <= highest_n_edits:
                    new_RAlist.append( read )
            self.RAlists[ i ] = new_RAlist

        return

    cpdef n_edits_sum ( self ):
        """
        """
        cdef:
            list n_edits_list
            object read
            int highest_n_edits

        n_edits_list = []

        for ralist in self.RAlists:
            for read in ralist:
                n_edits_list.append( read["n_edits"] )

        n_edits_list.sort()
        print ( n_edits_list )
        c = Counter( n_edits_list )
        print( c )

    cdef tuple __get_peak_REFSEQ ( self ):
        """Get the reference sequence within the peak region.

        """
        cdef:
            bytearray peak_refseq
            int i
            long prev_r                   #remember the previous filled right end
            long start
            long end
            long ind, ind_r
            object read
            bytearray read_refseq_ext
            bytearray read_refseq

        start = min( self.RAs_left, self.left )
        end = max( self.RAs_right, self.right )
        #print ("left",start,"right",end)
        peak_refseq_ext = bytearray( b'N' * ( end - start ) )

        # for treatment.
        peak_refseq_ext = self.__fill_refseq ( peak_refseq_ext, self.RAlists[0] )
        # and control if available.
        if self.RAlists[1]:
            peak_refseq_ext = self.__fill_refseq ( peak_refseq_ext, self.RAlists[1] )

        # trim
        peak_refseq = peak_refseq_ext[ self.left - start: self.right - start ]
        return ( bytes( peak_refseq ), bytes( peak_refseq_ext ) )

    cdef bytearray __fill_refseq ( self, bytearray seq, list ralist ):
        """Fill refseq sequence of whole peak with refseq sequences of
        each read in ralist.

        """
        cdef:
            long prev_r         # previous right position of last
                                # filled
            long ind, ind_r
            long start, end
            object read
            bytearray read_refseq

        start = min( self.RAs_left, self.left )

        prev_r = ralist[0]["lpos"]

        for i in range(  len( ralist ) ):
            read = ralist[ i ]
            if read[ "lpos" ] > prev_r:
                read = ralist[ i - 1 ]
                read_refseq = read.get_REFSEQ()
                ind = read["lpos"] - start
                ind_r = ind + read["rpos"] - read["lpos"]
                seq[ ind: ind_r  ] = read_refseq
                prev_r = read[ "rpos" ]
        # last
        read = ralist[ -1 ]
        read_refseq = read.get_REFSEQ()
        ind = read["lpos"] - start
        ind_r = ind + read["rpos"] - read["lpos"]
        seq[ ind: ind_r  ] = read_refseq
        return seq
                
    cpdef get_PosReadsInfo_ref_pos ( self, long ref_pos, bytes ref_nt ):
        """Generate a PosReadsInfo object for a given reference genome
        position.

        Return a PosReadsInfo object.

        """
        cdef:
            bytearray s
            bytearray bq
            object ra
            int i

        posreadsinfo_p = PosReadsInfo( ref_pos, ref_nt )
        #Treatment group
        for i in range( len( self.RAlists[ 0 ] ) ):
            ra = self.RAlists[ 0 ][ i ]
            if ra[ "lpos" ] <= ref_pos and ra[ "rpos" ] > ref_pos:
                ( s, bq ) = ra.get_variant_bq_by_ref_pos( ref_pos )
                if s == b'=':
                    posreadsinfo_p.add_T( i, ref_nt, bq[0] )
                elif len(s) < 2 and s != b'-':  # we will deal with SNV first
                    posreadsinfo_p.add_T( i, bytes(s[0:1]), bq[0] )

        #Control group
        for i in range( len( self.RAlists[ 1 ] ) ):
            ra = self.RAlists[ 1 ][ i ]
            if ra[ "lpos" ] <= ref_pos and ra[ "rpos" ] > ref_pos:
                ( s, bq ) = ra.get_variant_bq_by_ref_pos( ref_pos )
                if s == b'=':
                    posreadsinfo_p.add_C( i, ref_nt, bq[0] )
                elif len(s) < 2 and s != b'-':  # we will deal with SNV first
                    posreadsinfo_p.add_C( i, bytes(s[0:1]), bq[0] )
        return posreadsinfo_p

    cpdef bytearray get_FASTQ ( self ):
        """Get FASTQ file for all reads in RACollection.

        """
        cdef:
            object ra
            bytearray fastq_text
        
        fastq_text = bytearray(b"")

        for ra in self.RAlists[0]:
            fastq_text += ra.get_FASTQ()

        for ra in self.RAlists[1]:
            fastq_text += ra.get_FASTQ()
        
        return fastq_text

    cpdef list fermi_assemble( self, float fermiOverlapMinRatio ):
        """A wrapper function to call Fermi unitig building functions.
        """
        cdef:
            fml_opt_t *opt
            int c, n_seqs
            int * n_utg
            bseq1_t *seqs
            fml_utg_t *utg
            fml_utg_t p

            int unitig_k, merge_min_len
            bytes tmps
            bytes tmpq
            int ec_k = -1
            int64_t l
            bytes seq  #contains sequences of ALL reads, separated by '\x00';
            bytes qual #contains quality string of ALL reads, separated by '\x00';
            char * cseq
            char * cqual
            int i, j
            bytes tmpunitig
            bytes unitig                 #final unitig
            list unitig_list             # contain list of sequences in bytes format
            
        seq = b''
        qual = b''

        unitig_k=int(self.RAlists[0][0]["l"]*fermiOverlapMinRatio)
        merge_min_len=int(self.RAlists[0][0]["l"]*0.75)+1;

        n_seqs = len(self.RAlists[0]) + len(self.RAlists[1])

        # prepare seq and qual, note, we only extract SEQ according to the +
        # strand of reference sequence.
        seqs = <bseq1_t *> malloc( n_seqs * sizeof(bseq1_t) ) # we rely on fermi-lite to free this mem
        #print ("haha1")
        
        i = 0
        for ra in self.RAlists[0]:
            tmps = ra["SEQ"]
            tmpq = ra["QUAL"]
            l = len(tmps)

            cseq = <char *>malloc( (l+1)*sizeof(char))# we rely on fermi-lite to free this mem
            cqual = <char *>malloc( (l+1)*sizeof(char))# we rely on fermi-lite to free this mem
            for j in range(l):
                cseq[ j ] = tmps[ j ]
                cqual[ j ] = tmpq[ j ]
            cseq[ j ] = b'\x00'
            cqual[ j ]= b'\x00'

            seqs[ i ].seq = cseq
            seqs[ i ].qual = cqual
            seqs[ i ].l_seq = len(tmps)
            i += 1

        for ra in self.RAlists[1]:
            tmps = ra["SEQ"]
            tmpq = ra["QUAL"]
            l = len(tmps)

            cseq = <char *>malloc( (l+1)*sizeof(char))# we rely on fermi-lite to free this mem
            cqual = <char *>malloc( (l+1)*sizeof(char))# we rely on fermi-lite to free this mem
            for j in range(l):
                cseq[ j ] = tmps[ j ]
                cqual[ j ] = tmpq[ j ]
            cseq[ j ] = b'\x00'
            cqual[ j ]= b'\x00'

            seqs[ i ].seq = cseq
            seqs[ i ].qual = cqual
            seqs[ i ].l_seq = len(tmps)
            i += 1
        
        opt = <fml_opt_t *> PyMem_Malloc( sizeof(fml_opt_t) )
        opt.ec_k = -1
        opt.min_asm_ovlp = unitig_k
        opt.min_merge_len = merge_min_len
        #print ("haha2")

        n_utg = <int *> PyMem_Malloc( sizeof(int) )
        #print ("haha3")

        fml_opt_init(opt)
        #print ("haha3.1")
        utg = fml_assemble(opt, n_seqs, seqs, n_utg)
        #print ("haha4")        
        # get results
        unitig_list = []
        for i in range( n_utg[0] ):
            p = utg[ i ]
            if (p.len < 0):
                continue
            #unitig = b''
            #for j in range( p.len ):
            #    unitig += [b'A',b'C',b'G',b'T',b'N'][int(p.seq[j]) - 1]
            #unitig_list.append( unitig )
            unitig_list.append( p.seq )

        #print ("haha5")
        fml_utg_destroy(n_utg[0], utg)

        #PyMem_Free( seqs )
        #print ("haha6")
        PyMem_Free( opt )
        PyMem_Free( n_utg )
        #del seqs
        #del opt
        #del n_utg

        return unitig_list

    # cpdef list fermi_assemble( self, float fermiOverlapMinRatio ):
    #     """A wrapper function to call Fermi unitig building functions.
    #     """
    #     cdef:
    #         int unitig_k, merge_min_len
    #         bytes tmps
    #         bytes tmpq
    #         int ec_k = -1
    #         int64_t l
    #         mag_t *g
    #         magv_t p
    #         magopt_t *opt
    #         bytes seq  #contains sequences of ALL reads, separated by '\x00';
    #         bytes qual #contains quality string of ALL reads, separated by '\x00';
    #         char * cseq
    #         char * cqual
    #         int i
    #         kstring_t out
    #         bytes tmpunitig
    #         bytes unitig                 #final unitig
    #         list unitig_list             # contain list of sequences in bytes format
            
    #     seq = b''
    #     qual = b''

    #     unitig_k=int(self.RAlists[0][0]["l"]*fermiOverlapMinRatio)
    #     merge_min_len=int(self.RAlists[0][0]["l"]*0.75)+1;
        
    #     # prepare seq and qual, note, we only extract SEQ according to the +
    #     # strand of reference sequence.
    #     for ra in self.RAlists[0]:
    #         tmps = ra["SEQ"]
    #         tmpq = ra["QUAL"]
    #         seq += tmps + b'\x00'
    #         qual+= tmpq + b'\x00'

    #     for ra in self.RAlists[1]:
    #         tmps = ra["SEQ"]
    #         tmpq = ra["QUAL"]
    #         seq += tmps + b'\x00'
    #         qual+= tmpq + b'\x00'
            
    #     l = len( seq )
    #     cseq = seq
    #     cqual = qual

    #     # correct seq with qual
    #     fm6_api_correct(ec_k, l, cseq, cqual)
    #     # assemble unitigs
    #     g = fm6_api_unitig(unitig_k, l, cseq)
    #     opt = mag_init_opt()
    #     opt.flag |= 0x10 #MOG_F_CLEAN
    #     opt.min_merge_len = merge_min_len
    #     # clean
    #     mag_g_clean(g, opt)
    #     free(opt)

    #     # get results

    #     unitig_list = []
    #     for i in range( g.v.n ):
    #         p = g.v.a[ i ]
    #         if (p.len < 0):
    #             continue
    #         unitig = b''
    #         for j in range( p.len ):
    #             unitig += [b'A',b'C',b'G',b'T',b'N'][int(p.seq[j]) - 1]
    #         unitig_list.append( unitig )

    #     mag_g_destroy(g)
    #     return unitig_list

    cpdef tuple align_unitig_to_REFSEQ ( self, list unitig_list ):
        """Note: we use smith waterman, but we don't use linear gap
        penalty at this time.
        """
        cdef:
            bytes unitig
            seq_pair_t problem
            char * tmp
            bytes target
            bytes reference
            bytes target_aln
            bytes reference_aln
            list target_alns = []
            list reference_alns = []

        reference = copy(self.peak_refseq_ext+b'\x00')

        for unitig in unitig_list:
            target = copy(unitig + b'\x00')

            # we use swalign.c for local alignment (without affine gap
            # penalty). Will revise later.
            problem.a = target
            problem.alen = len( unitig )
            problem.b = reference
            problem.blen = len( self.peak_refseq_ext )
            results = smith_waterman( &problem )
            target_aln = results.seqs.a
            reference_aln = results.seqs.b
            # end of local alignment

            target_alns.append( target_aln )
            reference_alns.append( reference_aln )
            
        return ( target_alns, reference_alns )


    cpdef object remap_RAs_w_unitigs ( self, list unitig_list, tuple alns ):
        """unitig_list and tuple_alns are in the same order!

        return UnitigCollection object.

        """
        cdef:
            long start, end
            list target_alns, reference_alns
            list RAlists_T = []
            list RAlists_C = []
            object tmp_ra
            bytes tmp_ra_seq
            bytes tmp_unitig_seq
            bytes tmp_reference_seq
            bytes tmp_unitig_aln
            bytes tmp_reference_aln
            int i, j
            long left_padding_ref, right_padding_ref
            long left_padding_unitig, right_padding_unitig
            list ura_list = []
            object unitig_collection

        ( target_alns, reference_alns ) = alns

        start = min( self.left, self.RAs_left )
        end = max( self.right, self.RAs_right )
        
        for i in range( len(unitig_list) ):
            RAlists_T.append([])         # for each unitig, there is another list of RAs
            RAlists_C.append([])

        # assign RAs to unitigs
        for tmp_ra in self.RAlists[0]:
            tmp_ra_seq = tmp_ra["SEQ"]
            for i in range( len(unitig_list) ):
                unitig = unitig_list[ i ]
                if tmp_ra_seq in unitig:
                    RAlists_T[ i ].append( tmp_ra )
                    break
        for tmp_ra in self.RAlists[1]:
            tmp_ra_seq = tmp_ra["SEQ"]
            for i in range( len(unitig_list) ):
                unitig = unitig_list[ i ]
                if tmp_ra_seq in unitig:
                    RAlists_C[ i ].append( tmp_ra )
                    break

        # create UnitigCollection
        for i in range( len( unitig_list ) ):
            #b'---------------------------AAATAATTTTATGTCCTTCAGTACAAAAAGCAGTTTCAACTAAAACCCAGTAACAAGCTAGCAATTCCTTTTAAATGGTGCTACTTCAAGCTGCAGCCAGGTAGCTTTTTATTACAAAAAATCCCACAGGCAGCCACTAGGTGGCAGTAACAGGCTTTTGCCAGCGGCTCCAGTCAGCATGGCTTGACTGTGTGCTGCAGAAACTTCTTAAATCGTCTGTGTTTGGGACTCGTGGGGCCCCACAGGGCTTTACAAGGGCTTTTTAATTTCCAAAAACATAAAACAAAAAAA--------------'
            #b'GATATAAATAGGATGTTATGAGTTTTCAAATAATTTTATGTCCTTCAGTACAAAAAGCAGTTTCAACTAAAACCCAGTAACAAGCTAGCAATTCCTTTTAAATGGTGCTACTTCAAGCTGCAGCCAGGTAGCTTTTTATTACAAAAA-TCCCACAGGCAGCCACTAGGTGGCAGTAACAGGCTTTTGCCAGCGGCTCCAGTCAGCATGGCTTGACTGTGTGCTGCAGAAACTTCTTAAATCGTCTGTGTTTGGGACTCGTGGGGCCCCACAGGGCTTTACAAGGGCTTTTTAATTTCCAAAAACATAAAACAAAAAAAAATACAAATGTATT'
            tmp_unitig_aln = alns[ 0 ][ i ]
            tmp_reference_aln = alns[ 1 ][ i ]
            tmp_unitig_seq = tmp_unitig_aln.replace(b'-',b'')
            tmp_reference_seq = tmp_reference_aln.replace(b'-',b'')
            
            # print "tmp unitig aln:", tmp_unitig_aln
            # print "tmp refere aln:", tmp_reference_aln
            # print "peak refseqext:", self.peak_refseq_ext

            # find the position on self.peak_refseq_ext
            left_padding_ref = self.peak_refseq_ext.find( tmp_reference_seq ) # this number of nts should be skipped on refseq_ext from left
            right_padding_ref = len(self.peak_refseq_ext) - left_padding_ref - len(tmp_reference_seq) # this number of nts should be skipped on refseq_ext from right
            
            # print "padding ref", (left_padding_ref, right_padding_ref)

            #now, decide the lpos and rpos on reference of this unitig
            #first, trim left padding '-'
            left_padding_unitig = len(tmp_unitig_aln) - len(tmp_unitig_aln.lstrip(b'-'))
            right_padding_unitig = len(tmp_unitig_aln) - len(tmp_unitig_aln.rstrip(b'-'))

            tmp_lpos = start + left_padding_ref
            tmp_rpos = end - right_padding_ref

            for j in range( left_padding_unitig ):
                if tmp_reference_aln[ j ] != b'-':
                    tmp_lpos += 1
            for j in range( 1, right_padding_unitig + 1 ):
                if tmp_reference_aln[ -j ] != b'-':
                    tmp_rpos -= 1

            # print "unitig lpos and rpos:", ( tmp_lpos, tmp_rpos )


            tmp_unitig_aln = tmp_unitig_aln[ left_padding_unitig:(len(tmp_unitig_aln)-right_padding_unitig)]
            tmp_reference_aln = tmp_reference_aln[ left_padding_unitig:(len(tmp_reference_aln)-right_padding_unitig)]

            ura_list.append( UnitigRAs( self.chrom, tmp_lpos, tmp_rpos, tmp_unitig_aln, tmp_reference_aln, [RAlists_T[i], RAlists_C[i]] ) )

        return UnitigCollection( self.chrom, self.peak, ura_list )
                
