# Time-stamp: <2017-09-26 17:15:42 Tao Liu>

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

    cpdef bool is_indel ( self ):
        if self.v_mutation_type.find("Insertion") != -1 or self.v_mutation_type.find("Deletion") != -1:
            return True
        
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
        list l_chrom
        list l_p
        list l_Variants

    def __init__ ( self ):
        self.l_Variants = []
        self.l_chrom = []
        self.l_p = []

    cpdef add_Variant ( self, str chrom, long p, Variant v ):
        self.l_chrom.append( chrom )
        self.l_p.append( p )
        self.l_Variants.append( v )

    cpdef bool has_indel ( self ):
        cdef:
            int i
        for i in range( len(self.l_Variants) ):
            if self.l_Variants[ i ].is_indel():
                return True
        return False

    cpdef str toVCF ( self ):
        cdef:
            int i
            str res
        res = ""
        for i in range( len(self.l_Variants) ):
            res += "\t".join( ( self.l_chrom[ i ], str(self.l_p[ i ]+1), ".", self.l_Variants[ i ].toVCF() ) ) + "\n"
        return res
            
        
