# Time-stamp: <2017-11-01 13:13:13 Tao Liu>

"""Module for SAPPER PeakVariants class.

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
from SAPPER.Constants import *
from cpython cimport bool

cdef class Variant:
    cdef:
        long v_ref_pos
        str v_ref_allele
        str v_alt_allele
        int v_GQ
        str v_filter
        str v_type
        str v_mutation_type
        str v_top1allele
        str v_top2allele
        int v_DPT
        int v_DPC
        int v_DP1T
        int v_DP2T
        int v_DP1C
        int v_DP2C
        int v_PLUS1T
        int v_PLUS2T
        int v_MINUS1T
        int v_MINUS2T
        float v_deltaBIC
        float v_BIC_homo_major
        float v_BIC_homo_minor
        float v_BIC_heter_noAS
        float v_BIC_heter_AS
        float v_AR
        str v_GT
        int v_DP
        int v_PL_00
        int v_PL_01
        int v_PL_11
        
    def __init__ ( self, str ref_allele, str alt_allele, int GQ, str filter, str type, str mutation_type,
                   str top1allele, str top2allele, int DPT, int DPC, int DP1T, int DP2T, int DP1C, int DP2C,
                   int PLUS1T, int PLUS2T, int MINUS1T, int MINUS2T, 
                   float deltaBIC, float BIC_homo_major, float BIC_homo_minor, float BIC_heter_noAS, float BIC_heter_AS,
                   float AR, str GT, int DP, int PL_00, int PL_01, int PL_11):
        #self.v_ref_pos = ref_pos
        self.v_ref_allele = ref_allele
        self.v_alt_allele = alt_allele
        self.v_GQ = GQ
        self.v_filter = filter
        self.v_type = type
        self.v_mutation_type = mutation_type
        self.v_top1allele = top1allele
        self.v_top2allele = top2allele
        self.v_DPT = DPT
        self.v_DPC = DPC
        self.v_DP1T = DP1T
        self.v_DP2T = DP2T
        self.v_DP1C = DP1C
        self.v_DP2C = DP2C
        self.v_PLUS1T = PLUS1T
        self.v_PLUS2T = PLUS2T
        self.v_MINUS1T = MINUS1T
        self.v_MINUS2T = MINUS2T
        self.v_deltaBIC = deltaBIC
        self.v_BIC_homo_major = BIC_homo_major
        self.v_BIC_homo_minor = BIC_homo_minor
        self.v_BIC_heter_noAS = BIC_heter_noAS
        self.v_BIC_heter_AS = BIC_heter_AS
        self.v_AR = AR
        self.v_GT = GT
        self.v_DP = DP
        self.v_PL_00 = PL_00
        self.v_PL_01 = PL_01
        self.v_PL_11 = PL_11

    def __getstate__ ( self ):
        return (
            #self.v_ref_pos,
            self.v_ref_allele,
            self.v_alt_allele,
            self.v_GQ,
            self.v_filter,
            self.v_type,
            self.v_mutation_type,
            self.v_top1allele,
            self.v_top2allele,
            self.v_DPT,
            self.v_DPC,
            self.v_DP1T,
            self.v_DP2T,
            self.v_DP1C,
            self.v_DP2C,
            self.v_PLUS1T,
            self.v_PLUS2T,
            self.v_MINUS1T,
            self.v_MINUS2T,
            self.v_deltaBIC,
            self.v_BIC_homo_major,
            self.v_BIC_homo_minor,
            self.v_BIC_heter_noAS,
            self.v_BIC_heter_AS,
            self.v_AR,
            self.v_GT,
            self.v_DP,
            self.v_PL_00,
            self.v_PL_01,
            self.v_PL_11 )

    def __setstate__ ( self, state ):
        ( #self.v_ref_pos,
          self.v_ref_allele,
          self.v_alt_allele,
          self.v_GQ,
          self.v_filter,
          self.v_type,
          self.v_mutation_type,
          self.v_top1allele,
          self.v_top2allele,
          self.v_DPT,
          self.v_DPC,
          self.v_DP1T,
          self.v_DP2T,
          self.v_DP1C,
          self.v_DP2C,
          self.v_PLUS1T,
          self.v_PLUS2T,
          self.v_MINUS1T,
          self.v_MINUS2T,
          self.v_deltaBIC,
          self.v_BIC_homo_major,
          self.v_BIC_homo_minor,
          self.v_BIC_heter_noAS,
          self.v_BIC_heter_AS,
          self.v_AR,
          self.v_GT,
          self.v_DP,
          self.v_PL_00,
          self.v_PL_01,
          self.v_PL_11 ) = state

    cpdef bool is_indel ( self ):
        if self.v_mutation_type.find("Insertion") != -1 or self.v_mutation_type.find("Deletion") != -1:
            return True
        else:
            return False

    cpdef bool only_del ( self ):
        if self.v_mutation_type == "Deletion":
            return True
        else:
            return False

    def __getitem__ ( self, keyname ):
        if keyname == "ref_allele":
            return self.v_ref_allele
        elif keyname == "alt_allele":
            return self.v_alt_allele
        elif keyname == "type":
            return self.type
        elif keyname == "mutation_type":
            return self.mutation_type
        else:
            raise Exception("keyname is not accessible:", keyname)

    def __setitem__ ( self, keyname, v ):
        if keyname == "ref_allele":
            self.v_ref_allele = v
        elif keyname == "alt_allele":
            self.v_alt_allele = v
        elif keyname == "type":
            self.type = v
        elif keyname == "mutation_type":
            self.mutation_type = v
        else:
            raise Exception("keyname is not accessible:", keyname)        
        
    cpdef bool is_refer_biased_01 ( self, float ar=0.85 ):
        if self.v_AR >= ar and self.v_ref_allele == self.v_top1allele:
            return True
        else:
            return False
        
    cpdef str toVCF ( self ):
        return "\t".join( ( self.v_ref_allele, self.v_alt_allele, "%d" % self.v_GQ, self.v_filter,
                            "M=%s;MT=%s;DPT=%d;DPC=%d;DP1T=%d%s;DP2T=%d%s;DP1C=%d%s;DP2C=%d%s;SB=%d,%d,%d,%d;DBIC=%.2f;BICHOMOMAJOR=%.2f;BICHOMOMINOR=%.2f;BICHETERNOAS=%.2f;BICHETERAS=%.2f;AR=%.2f" % \
                                (self.v_type, self.v_mutation_type, self.v_DPT, self.v_DPC, self.v_DP1T, self.v_top1allele,
                                 self.v_DP2T, self.v_top2allele, self.v_DP1C, self.v_top1allele, self.v_DP2C, self.v_top2allele,
                                 self.v_PLUS1T, self.v_PLUS2T, self.v_MINUS1T, self.v_MINUS2T,
                                 self.v_deltaBIC,
                                 self.v_BIC_homo_major, self.v_BIC_homo_minor, self.v_BIC_heter_noAS,self.v_BIC_heter_AS,
                                 self.v_AR
                                 ),
                            "GT:DP:GQ:PL",
                            "%s:%d:%d:%d,%d,%d" % (self.v_GT, self.v_DP, self.v_GQ, self.v_PL_00, self.v_PL_01, self.v_PL_11)
                            ) )

cdef class PeakVariants:
    cdef:
        str chrom
        dict d_Variants

    def __init__ ( self, str chrom ):
        self.chrom = chrom
        self.d_Variants = {}

    def __getstate__ ( self ):
        return ( self.d_Variants, self.chrom )

    def __setstate__ ( self, state ):
        ( self.d_Variants, self.chrom ) = state

    cpdef int n_variants ( self ):
        return len(self.d_Variants)
            
    cpdef add_variant ( self, long p, Variant v ):
        self.d_Variants[ p ] = v

    cpdef bool has_indel ( self ):
        cdef:
            long p
        for p in sorted( self.d_Variants.keys() ):
            if self.d_Variants[ p ].is_indel():
                return True
        return False

    cpdef bool has_refer_biased_01 ( self ):
        cdef:
            long p 
        for p in sorted( self.d_Variants.keys() ):
            if self.d_Variants[ p ].is_refer_biased_01():
                return True
        return False

    cpdef list get_refer_biased_01s ( self ):
        cdef:
            list ret_poss = []
            long p
        for p in sorted( self.d_Variants.keys() ):
            if self.d_Variants[ p ].is_refer_biased_01():
                ret_poss.append( p )
        return ret_poss

    cpdef remove_variant ( self, long p ):
        assert p in self.d_Variants
        self.d_Variants.pop( p )
    
    cpdef replace_variant ( self, long p, Variant v ):
        assert p in self.d_Variants
        self.d_Variants[ p ] = v


    cpdef merge_continuous_dels ( self ):
        cdef:
            long p0, p
        p0 = -1
        for p in sorted( self.d_Variants.keys() ):
            if p == p0+1 and self.d_Variants[ p ].only_del() and self.d_Variants[ p0 ].only_del() :
                # we keep p0, remove p, and add p's ref_allele to p0, keep other information as in p0
                self.d_Variants[ p0 ]["ref_allele"] += self.d_Variants[ p ]["ref_allele"]
                self.d_Variants.pop ( p )
            p0 = p
        return
                
    cpdef str toVCF ( self ):
        cdef:
            long p
            str res
        res = ""
        for p in sorted( self.d_Variants.keys() ):
            res += "\t".join( ( self.chrom, str(p+1), ".", self.d_Variants[ p ].toVCF() ) ) + "\n"
        return res
            
        
