from core.MPIUtil cimport MPIUtil
from includes.file cimport FILE

cdef extern from "core/Util.h" nogil:
    MPIUtil* mpiUtil
    FILE* globalLog
    FILE* nullLog
    void initSystem(int argc, char** argv)
