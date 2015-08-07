cdef extern from "core/MPIUtil.h" nogil:
    cdef cppclass MPIUtil:
        MPIUtil(int argc, char** argv)
