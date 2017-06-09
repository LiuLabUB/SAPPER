# Time-stamp: <2017-06-09 16:46:59 Tao Liu>

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
from operator import itemgetter

from SAPPER.Constants import *
from SAPPER.ReadAlignment import ReadAlignment
from SAPPER.PosReadsInfo import PosReadsInfo
from SAPPER.PeakIO import PeakIO

from cpython cimport bool

import numpy as np
cimport numpy as np
from numpy cimport uint32_t, uint64_t, int32_t, int64_t

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

cdef extern from "kstring.h":
    ctypedef struct kstring_t:
        uint32_t l, m
        char *s

cdef extern from "priv.h":
    int ksa_bwt64(unsigned char *T, int64_t n, int k)

cdef extern from "mag.h":
    ctypedef struct magopt_t:
        int flag, max_arc, n_iter, min_ovlp, min_elen, min_ensr, min_insr, max_bdist, max_bvtx, min_merge_len
        float min_dratio0, min_dratio1
        float max_bcov, max_bfrac
    ctypedef struct ku128_t:
        uint64_t x, y
    ctypedef struct ku64_v:
        size_t n, m
        uint64_t *a
    ctypedef struct ku128_v:
        size_t n, m
        ku128_t *a
    ctypedef struct magv_t:
        int len, nsr            #length; number supporting reads
        uint32_t max_len        # allocated seq/cov size
        uint64_t k[2]           # bi-interval
        ku128_v nei[2]          # neighbors
        char *seq, *cov         # sequence and coverage
        void *ptr               # additional information
    ctypedef struct magv_v:
        size_t n, m
        magv_t *a
    ctypedef struct mag_t:
        magv_v v
        float rdist             # read distance
        int min_ovlp            # minimum overlap seen from the graph
        void *h

    magopt_t *mag_init_opt()
    void mag_g_clean(mag_t *g, const magopt_t *opt)
    void mag_g_print(const mag_t *g)
    void mag_g_destroy(mag_t *g)
    void mag_v_write(const magv_t *p, kstring_t *out);
    
cdef extern from "fermi.h":
    int fm6_api_correct(int kmer, int64_t l, char *_seq, char *_qual)
    mag_t *fm6_api_unitig(int min_match, int64_t l, char *seq)


    
# ------------------------------------
# constants
# ------------------------------------
__version__ = "Parser $Revision$"
__author__ = "Tao Liu <tliu4@buffalo.edu>"
__doc__ = "All Parser classes"

__DNACOMPLEMENT__ = b'\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQRSAUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff' # A trans table to convert A to T, C to G, G to C, and T to A.

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
        self.sorted = True
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
        else:
            raise KeyError("Unavailable key:", keyname)

    def __getstate__ ( self ):
        #return {"chrom":self.chrom, "peak":self.peak, "RAlists":self.RAlists,
        #        "left":self.left, "right":self.right, "length": self.length,
        #        "RAs_left":self.RAs_left, "RAs_right":self.RAs_right}
        return (self.chrom, self.peak, self.RAlists, self.left, self.right, self.length, self.RAs_left, self.RAs_right)
        
    def __setstate__ ( self, state ):
        (self.chrom, self.peak, self.RAlists, self.left, self.right, self.length, self.RAs_left, self.RAs_right) = state
        
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

    cpdef bytearray get_peak_REFSEQ ( self ):
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
            bytearray read_refseq

        start = min( self.RAs_left, self.left )
        end = max( self.RAs_right, self.right )
        peak_refseq = bytearray( b'N' * ( end - start ) )

        # for treatment.
        peak_refseq = self.__fill_refseq ( peak_refseq, self.RAlists[0] )
        # and control if available.
        if self.RAlists[1]:
            peak_refseq = self.__fill_refseq ( peak_refseq, self.RAlists[1] )

        # trim
        peak_refseq = peak_refseq[ self.left - start: self.right - start ]
        return peak_refseq

    cdef bytearray __fill_refseq ( self, seq, ralist ):
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

        start = self.RAs_left

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
            int unitig_k, merge_min_len
            bytes tmps
            bytes tmpq
            int ec_k = -1
            int64_t l
            mag_t *g
            magv_t p
            magopt_t *opt
            bytes seq  #contains sequences of ALL reads, separated by '\x00';
            bytes qual #contains quality string of ALL reads, separated by '\x00';
            char * cseq
            char * cqual
            int i
            kstring_t out
            bytes tmpunitig
            bytes unitig                 #final unitig
            list unitig_list
            
        seq = b''
        qual = b''

        unitig_k=int(self.RAlists[0][0]["l"]*fermiOverlapMinRatio)
        merge_min_len=int(self.RAlists[0][0]["l"]*0.75)+1;
        
        # prepare seq and qual
        for ra in self.RAlists[0]:
            ( tmps, tmpq ) =  ra.get_SEQ_QUAL()
            seq += tmps + b'\x00'
            qual+= tmpq + b'\x00'

        for ra in self.RAlists[1]:
            ( tmps, tmpq ) =  ra.get_SEQ_QUAL()
            seq += tmps + b'\x00'
            qual+= tmpq + b'\x00'
            
        l = len( seq )
        cseq = seq
        cqual = qual

        # correct seq with qual
        fm6_api_correct(ec_k, l, cseq, cqual)
        # assemble unitigs
        g = fm6_api_unitig(unitig_k, l, cseq)
        opt = mag_init_opt()
        opt.flag |= 0x10 #MOG_F_CLEAN
        opt.min_merge_len = merge_min_len
        # clean
        mag_g_clean(g, opt)
        free(opt)

        # get results

        unitig_list = []
        for i in range( g.v.n ):
            p = g.v.a[ i ]
            if (p.len < 0):
                continue
            unitig = b''
            for j in range( p.len ):
                unitig += [b'A',b'C',b'G',b'T',b'N'][int(p.seq[j]) - 1]
            unitig_list.append( unitig )

        mag_g_destroy(g)
        return unitig_list

    cpdef align_unitig_to_REFSEQ ( self, list unitig_list, bytes refseq ):
        """
        """
        cdef:
            list unitig_refstart_list
            list unitig_refend_list
            
        return


    cpdef filter_RAs_w_unitigs ( self, list unitig_list ):
        """
        """
        cdef:
            int i
            bytes tmps
            bytes unitig
            list new_RAlist_T, new_RAlist_C

        new_RAlist_T = []
        new_RAlist_C = []
        for ra in self.RAlists[0]:
            tmps =  bytes(ra.get_SEQ())
            for unitig in unitig_list:
                if tmps in unitig:
                    new_RAlist_T.append( ra )
                    break
                else:
                    tmps = tmps[::-1].translate( __DNACOMPLEMENT__ )
                    if tmps in unitig:
                        new_RAlist_T.append( ra )
                        break

        for ra in self.RAlists[1]:
            tmps =  bytes(ra.get_SEQ())
            for unitig in unitig_list:
                if tmps in unitig:
                    new_RAlist_C.append( ra )
                    break
                else:
                    tmps = tmps[::-1].translate( __DNACOMPLEMENT__ )
                    if tmps in unitig:
                        new_RAlist_C.append( ra )
                        break

        self.RAlists = [new_RAlist_T, new_RAlist_C]
        self.RAs_left = new_RAlist_T[0]["lpos"] # initial assignment of RAs_left
        self.RAs_right = new_RAlist_T[-1]["rpos"] # initial assignment of RAs_right
        # check RAs_left and RAs_right
        for ra in new_RAlist_T:
            if ra[ "lpos" ] < self.RAs_left:
                self.RAs_left = ra[ "lpos" ]
            if ra[ "rpos" ] > self.RAs_right:
                self.RAs_right = ra[ "rpos" ]

        for ra in new_RAlist_C:
            if ra[ "lpos" ] < self.RAs_left:
                self.RAs_left = ra[ "lpos" ]
            if ra[ "rpos" ] > self.RAs_right:
                self.RAs_right = ra[ "rpos" ]
