from fluid.FluidParams cimport FluidSolverParams

cdef extern from "electronic/ElecVars.h" nogil:
    cdef cppclass ElecVars:
        FluidSolverParams fluidParams
        int lcaoIter
        double lcaoTol
        double subspaceRotationFactor
