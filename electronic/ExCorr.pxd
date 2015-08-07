from libcpp cimport bool

from includes.memory cimport shared_ptr
from electronic.Everything cimport Everything

cdef extern from "electronic/ExCorr.h" nogil:
    cdef enum ExCorrType:
        ExCorrLDA_PZ
        ExCorrLDA_PW
        ExCorrLDA_PW_prec
        ExCorrLDA_VWN
        ExCorrLDA_Teter
        ExCorrGGA_PBE
        ExCorrGGA_PBEsol
        ExCorrGGA_PW91
        ExCorrMGGA_TPSS
        ExCorrMGGA_revTPSS
        ExCorrORB_GLLBsc
        ExCorrPOT_LB94
        ExCorrHYB_PBE0
        ExCorrHYB_HSE06
        ExCorrHYB_HSE12
        ExCorrHYB_HSE12s
        ExCorrHF

    cdef enum KineticType:
        KineticNone
        KineticTF
        KineticVW
        KineticPW91

    cdef struct IncludeTXC:
        bool T
        bool X
        bool C

        IncludeTXC(bool T=False, bool X=True, bool C=True)

    cdef cppclass ExCorr:
        ExCorr(ExCorrType exCorrType = ExCorrGGA_PBE,
               KineticType kineticType = KineticNone)
        ExCorrType exCorrType
        KineticType kineticType
