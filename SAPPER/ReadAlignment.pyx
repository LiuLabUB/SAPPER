# Time-stamp: <2017-06-06 11:39:21 Tao Liu>

"""Module for SAPPER ReadAlignment class

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
import gzip
import io
from cpython cimport bool

from SAPPER.Constants import *

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

__BAMDNACODE__ = b"=ACMGRSVTWYHKDBN"
__CIGARCODE__ = "MIDNSHP=X"

# -- CIGAR CODE --
#OP BAM  Description
#M  0    alignment match (can be a sequence match or mismatch) insertion to the reference
#I  1    insertion to the reference
#D  2    deletion from the reference
#N  3    skipped region from the reference
#S  4    soft clipping (clipped sequences present in SEQ)
#H  5    hard clipping (clipped sequences NOT present in SEQ)
#P  6    padding (silent deletion from padded reference)
#=  7    sequence match
#X  8    sequence mismatch
# -- -- -- -- -- --

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------
cdef class ReadAlignment:
    cdef:
        bytes chrom
        int lpos
        int rpos
        int strand
        bytes binaryseq
        bytes binaryqual
        tuple cigar             # each item contains op_l|op
        bytes MD
        int n_edits             # number of edits; higher the number,
                                # more differences with reference

    def __init__ ( self, 
                   bytes chrom, int lpos, int rpos,
                   int strand,
                   bytes binaryseq, 
                   bytes binaryqual,
                   tuple cigar,
                   bytes MD ):
        self.chrom = chrom
        self.lpos = lpos
        self.rpos = rpos
        self.strand = strand
        self.binaryseq = binaryseq
        self.binaryqual = binaryqual
        self.cigar = cigar
        self.MD = MD
        self.n_edits = self.get_n_edits()

    cdef int get_n_edits( self ):
        """The number is from self.cigar and self.MD.

        """
        cdef:
            int n_edits
            int i, cigar_op, cigar_op_l
            char c
        
        n_edits = 0
        for i in self.cigar:    # only count insertion or softclip
            cigar_op = i & 15
            cigar_op_l = i >> 4
            if cigar_op in [ 1, 4 ]:    # count Insertion or Softclip
                n_edits += cigar_op_l
        
        for c in self.MD:
            if (c > 64 and c < 91): # either deletion in query or mismatch
                n_edits += 1
        return n_edits

    def __getitem__ ( self, keyname ):
        if keyname == "chrom":
            return self.chrom
        elif keyname == "lpos":
            return self.lpos
        elif keyname == "rpos":
            return self.rpos
        elif keyname == "strand":
            return self.strand
        elif keyname == "binaryseq":
            return self.binaryseq
        elif keyname == "binaryqual":
            return self.binaryqual
        elif keyname == "cigar":
            return self.cigar
        elif keyname == "MD":
            return self.MD
        elif keyname == "n_edits":
            return self.n_edits
        else:
            raise KeyError("No such key", keyname)

    def __setitem__ ( self, keyname, value ):
        if keyname == "chrom":
            self.chrom = value
        elif keyname == "lpos":
            self.lpos = value
        elif keyname == "rpos":
            self.rpos = value
        elif keyname == "strand":
            self.strand = value
        # do not allow change to binaryseq
        #elif keyname == "binaryseq":
        #    self.binaryseq = value
        # allow change to bq, in case of bq calibration
        elif keyname == "binaryqual":
            self.binaryqual = value
        elif keyname == "cigar":
            self.cigar = value
        elif keyname == "MD":
            self.MD = value
        # do not allow change to n_edits
        #elif keyname == "n_edits":
        #    self.n_edits = value
        else:
            raise KeyError("No such key", keyname)

    def __getstate__ ( self ):
        return ( self.chrom, self.lpos, self.rpos, self.strand, self.binaryseq, self.binaryqual, self.cigar, self.MD, self.n_edits )

    def __setstate__ ( self, state ):
        ( self.chrom, self.lpos, self.rpos, self.strand, self.binaryseq, self.binaryqual, self.cigar, self.MD, self.n_edits ) = state


    cpdef bytearray get_SEQ ( self ):
        """Convert binary seq to ascii seq.

        Rule: for each byte, 1st base in the highest 4bit; 2nd in the lowest 4bit. "=ACMGRSVTWYHKDBN" -> [0,15]

        Note: In BAM, if a sequence is mapped to reverse strand, the
        reverse complement seq is written in SEQ field. So the return
        value of this function will not be the original one if the
        read is mapped to - strand.
        """
        cdef:
            char c
            bytearray seq
            bytes s
        seq = bytearray(b"")
        for c in self.binaryseq:
            # high
            seq.append( __BAMDNACODE__[c >> 4 & 15] )
            # low
            seq.append( __BAMDNACODE__[c & 15] )
        if seq[-1] == b"=":
            # trim the last '=' if it exists
            seq = seq[:-1]
        return seq

    cpdef bytearray get_REFSEQ ( self ):
        """Fetch reference sequence, using self.MD and self.cigar
        """
        cdef:
            char c
            bytearray seq, refseq
            int i, cigar_op, cigar_op_l
            bytearray MD_op
            int ind
            bool flag_del       # flag for deletion event in query

        seq = self.get_SEQ()    # we start with read seq then make modifications

        # 2-step proces
        # First step: use CIGAR to edit SEQ to remove S (softclip) and I (insert)
        # __CIGARCODE__ = "MIDNSHP=X"
        # let ind be the index in SEQ
        ind = 0
        for i in self.cigar:
            cigar_op = i & 15
            cigar_op_l = i >> 4
            if cigar_op in [2, 5, 6]:     # do nothing for Deletion (we will
                                          # put sequence back in step 2),
                                          # Hardclip and Padding
                pass
            elif cigar_op in [0, 7, 8]:   # M = X alignment match (match or
                                          # mismatch)
                # do nothing and move ind
                ind += cigar_op_l
            elif cigar_op in [ 1, 4 ]:    # Remove for Insertion or Softclip
                seq[ ind : ind + cigar_op_l ] = b''

        # now the seq should be at the same length as rpos-lpos 

        # Second step: use MD string to edit SEQ to put back 'deleted
        # seqs' and modify mismatches

        # let ind be the index in SEQ again, from 0
        ind = 0
        MD_op = bytearray(b'')
        flag_del = False
        for c in self.MD:
            if c < 58 and c > 47:
                # means Match
                flag_del = False
                MD_op.append(c)
            elif (c > 64 and c < 91) and not flag_del:
                # An alphabet means Mismatch, Note, if MD is made
                # right, a mismatch should only be 1 letter surrounded
                # by digits.
                ind += int(MD_op)
                seq[ ind ] = c
                ind += 1
                # reset MD_op
                MD_op = bytearray(b'')
            elif (c > 64 and c < 91) and flag_del:
                seq[ ind:ind ] = [c,]
                ind += 1
            elif c == 94:
                # means Deletion in query. Now, insert a sequnce into
                # SEQ
                flag_del = True
                ind += int(MD_op)
                # reset MD_op
                MD_op = bytearray(b'')
            else:
                raise Exception("Don't understand this operator in MD: %c" % c)
            #print( seq.decode() )

        return seq
    
    cpdef get_base_by_ref_pos ( self, long ref_pos ):
        """Get base by ref position.

        """
        cdef:
           int relative_pos, p
        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")
        relative_pos = ref_pos - self.lpos
        p = self.relative_ref_pos_to_relative_query_pos( relative_pos ) 

        if p == -1:             # located in a region deleted in query
            return None
        else:
            return __BAMDNACODE__[ (self.binaryseq[p//2] >> ((1-p%2)*4) ) & 15 ]

    cpdef get_bq_by_ref_pos ( self, long ref_pos ):
        """Get base quality by ref position. Base quality is in Phred scale.

        Returned value is the raw Phred-scaled base quality.

        """
        cdef:
           int relative_pos, p
        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")
        relative_pos = ref_pos - self.lpos
        p = self.relative_ref_pos_to_relative_query_pos( relative_pos ) 

        if p == -1:             # located in a region deleted in query
            return None
        else:
            return self.binaryqual[p]

    cpdef tuple get_base_bq_by_ref_pos ( self, long ref_pos ):
        """Get base and base quality by ref position. Base quality is in Phred scale.

        Returned bq is the raw Phred-scaled base quality.
        """
        cdef:
           int relative_pos, p
        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")
        relative_pos = ref_pos - self.lpos
        p = self.relative_ref_pos_to_relative_query_pos( relative_pos ) 

        if p == -1:             # located in a region deleted in query
            return None
        else:
            return ( __BAMDNACODE__[ (self.binaryseq[p//2] >> ((1-p%2)*4) ) & 15 ], self.binaryqual[p] )

    cpdef tuple get_variant_bq_by_ref_pos ( self, long ref_pos ):
        """Get any variants (different with reference) and base quality by ref position. 

        variants will be 

        1) =, if identical

        2) A/T/C/G, if SNV

        3) -, if the reference base is deleted, in this case, bq will
        be the highest possible bq, which is 93.

        4) ^<A/T/C/G>+, if there is an insertion at the location

        Base quality is the raw Phred-scaled base quality.

        """
        cdef:
           int res, p, op, op_l
           bytearray refseq
           bytes p_refseq, p_seq
           bytearray seq_array
           bytearray bq_array

        assert self.lpos <= ref_pos and self.rpos > ref_pos, Exception("Given position out of alignment location")

        res = ref_pos - self.lpos          # residue
        p = 0
        
        refseq = self.get_REFSEQ()
        p_refseq =  refseq[ res ]
        # -- CIGAR CODE --
        #OP BAM  Description
        #M  0    alignment match (can be a sequence match or mismatch) insertion to the reference
        #I  1    insertion to the reference
        #D  2    deletion from the reference
        #N  3    skipped region from the reference
        #S  4    soft clipping (clipped sequences present in SEQ)
        #H  5    hard clipping (clipped sequences NOT present in SEQ)
        #P  6    padding (silent deletion from padded reference)
        #=  7    sequence match
        #X  8    sequence mismatch

        seq_array = bytearray( b'' )
        bq_array = bytearray( b'' )

        #if ref_pos == 17255842 and self.lpos == 17255842:
        #    print( self.lpos, self.rpos, refseq, p_refseq )

        for i in self.cigar:
            op = i & 15
            op_l = i >> 4
            if op in [0, 7, 8]:         # M = X alignment match (match or mismatch)
                if res < op_l:
                    p += res
                    # find the position, now get the ref
                    seq_array.append( __BAMDNACODE__[ (self.binaryseq[ p//2 ] >> ( (1-p%2)*4 ) ) & 15 ] )
                    bq_array.append( self.binaryqual[ p ] )
                    if seq_array == p_refseq:
                        seq_array = bytearray( b'=' )
                        
                    #if ref_pos == 17255842 and self.lpos == 17255842:
                    #    print( p_refseq, seq_array, "mark" )

                    break
                else:
                    # go to the next cigar code
                    p += op_l
                    res -= op_l
            elif op in [ 2, 3 ]: # D N
                if res < op_l:
                    # find the position, however ...
                    # position located in a region in reference that not exists in query
                    seq_array.append( b'-' )
                    bq_array.append( 93 )

                    #if ref_pos == 17255842 and self.lpos == 17255842:
                    #    print( "mark2" )


                    break
                else:
                    # go to the next cigar code
                    res -= op_l
            elif op == 1 :      # Insertion
                if res == 0:    # no residue left, so return a chunk of inserted sequence
                    # first, add the insertion point
                    seq_array = bytearray( b'~' )
                    bq_array.append( self.binaryqual[ p ] )
                    # then add the inserted seq
                    for i in range( op_l ):
                        p += 1
                        seq_array.append( __BAMDNACODE__[ (self.binaryseq[ p//2 ] >> ( (1-p%2)*4 ) ) & 15 ]  ) 
                        bq_array.append( self.binaryqual[ p ] )

                    #if ref_pos == 17255842 and self.lpos == 17255842:
                    #    print( "mark3" )
                    
                    break
                else:
                    p += op_l
            elif op == 4 :      # Softclip. If it's Softclip, we'd better not return the extra seq
                p += op_l

        return ( seq_array, bq_array )
        # last position ?
        #raise Exception("Not expected to see this")

    cdef int relative_ref_pos_to_relative_query_pos ( self, long relative_ref_pos ):
        """Convert relative pos on ref to pos on query.
        """
        cdef:
            int p, res, op, op_l
        p = 0
        res = relative_ref_pos
        
        for i in self.cigar:
            op = i & 15
            op_l = i >> 4
            if op in [0, 7, 8]:         # M = X alignment match (match or mismatch)
                if res < op_l:
                    p += res
                    return p
                else:
                    p += op_l
                    res -= op_l
            elif op in [ 2, 3 ]: # D N
                if res < op_l:
                    # position located in a region in reference that not exists in query
                    return -1
                else:
                    res -= op_l
            elif op in [ 1, 4 ]:       # I
                p += op_l
        return p


### End ###

