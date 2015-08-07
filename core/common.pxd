#The corresponding file from jdftx lives under electronic but it has core stuff

cdef extern from "electronic/common.h" nogil:
    cdef cppclass matrix
    cdef cppclass diagMatrix
    cdef cppclass ColumnBundle
    cdef cppclass QuantumNumber
    cdef cppclass Control
    cdef cppclass GridInfo
    cdef cppclass Basis
    cdef cppclass SpeciesInfo
    cdef cppclass IonInfo
    cdef cppclass Symmetries
    cdef cppclass ExactExchange
    cdef cppclass ExCorr
    cdef cppclass ElecInfo
    cdef cppclass ElecVars
    cdef cppclass Energies
    cdef cppclass Everything

    cdef struct ElecGradient
    cdef struct IonicGradient
    cdef cppclass InverseKohnSham
    cdef cppclass VanDerWaals
    cdef cppclass Vibrations
