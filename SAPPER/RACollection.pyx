# Time-stamp: <2017-06-06 16:40:02 Tao Liu>

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
from numpy cimport uint32_t, uint64_t, int32_t

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

    cpdef get_peak_REFSEQ ( self ):
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

        start = self.RAs_left
        end = self.RAs_right
        peak_refseq = bytearray( b'N' * ( end - start ) )

        # prev_r = start

        # # first
        # read = self.RAlists[ 0 ][ 0 ]
        # read_refseq = read.get_REFSEQ()
        # ind = read["lpos"] - start
        # ind_r = ind + read["rpos"] - read["lpos"]
        # peak_refseq[ ind: ind_r  ] = read_refseq
        # prev_r = read[ "rpos" ]

        # # from 2nd read, skip if overlap
        # for i in range( 1, len( self.RAlists[ 0 ] ) ):
        #     read = self.RAlists[ 0 ][ i ]
        #     if read[ "lpos" ] > prev_r:
        #         read = self.RAlists[ 0 ][ i - 1 ]
        #         read_refseq = read.get_REFSEQ()
        #         ind = read["lpos"] - start
        #         ind_r = ind + read["rpos"] - read["lpos"]
        #         peak_refseq[ ind: ind_r  ] = read_refseq
        #         prev_r = read[ "rpos" ]
        # # last
        # read = self.RAlists[ 0 ][ -1 ]
        # read_refseq = read.get_REFSEQ()
        # ind = read["lpos"] - start
        # ind_r = ind + read["rpos"] - read["lpos"]
        # peak_refseq[ ind: ind_r  ] = read_refseq

        
        # for treatment.
        peak_refseq = self.__fill_refseq ( peak_refseq, self.RAlists[0] )
        # and control if available.
        if self.RAlists[1]:
            peak_refseq = self.__fill_refseq ( peak_refseq, self.RAlists[1] )

        # trim
        peak_refseq = peak_refseq[ self.left - self.RAs_left: self.right - self.RAs_left ]
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
