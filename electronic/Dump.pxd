from libcpp.set cimport set
from libcpp.pair cimport pair

cdef extern from "electronic/Dump.h" nogil:
    enum DumpFrequency:
        DumpFreq_End

    enum DumpVariable:
        DumpAll
        DumpNone

    cdef cppclass Dump(set[pair[DumpFrequency,DumpVariable]]):
        pass
