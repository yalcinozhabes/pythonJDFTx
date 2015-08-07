from libcpp cimport bool

from includes.vector3 cimport vector3
from includes.string cimport string

cdef extern from "electronic/Control.h" nogil:
    enum ElecEigenAlgo:
        ElecEigenCG
        ElecEigenDavidson

    enum BasisKdep:
        BasisKpointDep
        BasisKpointIndep

    cdef cppclass Control:
        bool fixed_H
        bool fixOccupied
        double occupiedThreshold
        bool cacheProjectors
        double davidsonBandRatio

        ElecEigenAlgo elecEigenAlgo
        BasisKdep basisKdep
        double Ecut, EcutRho

        bool dragWavefunctions
        vector3[double] lattMoveScale

        int fluidGummel_nIterations
        double fluidGummel_Atol

        double overlapConditionThreshold
        int overlapCheckInterval

        bool shouldPrintEigsFillings
        bool shouldPrintEcomponents
        bool shouldPrintMuSearch
        bool shouldPrintKpointsBasis

        bool invertKS
        bool invertKS_nonlocal
        double invertKS_sigma
        string invertKS_chiGuessFilename

        bool scf
        bool convergeEmptyStates
        bool dumpOnly

        Control()
