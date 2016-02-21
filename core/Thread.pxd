cdef extern from "core/Thread.h" nogil:
    int nProcsAvailable
    void resumeOperatorThreading()
