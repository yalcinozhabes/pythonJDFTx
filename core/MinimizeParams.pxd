from libcpp cimport bool
from includes.file cimport *

cdef extern from "core/MinimizeParams.h" nogil:
    cdef struct MinimizeParams:
        int nIterations
        int nDim
        int history
        FILE* fpLog
        const char* linePrefix
        const char* energyLabel
        const char* energyFormat
        double knormThreshold
        double energyDiffThreshold
        int nEnergyDiff

        double alphaTstart
        double alphaTmin
        bool updateTestStepSize

        double alphaTreduceFactor
        double alphaTincreaseFactor
        int nAlphaAdjustMax

        double wolfeEnergy
        double wolfeGradient

        bool fdTest
