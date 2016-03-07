cimport mpi4py.libmpi as mpi

from libcpp cimport bool
from libcpp.vector cimport vector

from core.common cimport *
from electronic.Everything cimport Everything
from includes.vector3 cimport vector3

cimport numpy as cnp

cdef extern from "electronic/ElecInfo.h" nogil:
    cdef cppclass QuantumNumber:
        vector3[double] k
        int spin
        double weight

        QuantumNumber()
        int index()

    cdef inline vector3[double] getCoord(QuantumNumber& qnum)

cdef void addQnum(Everything& e, cnp.ndarray kpt, double weight)
cdef void clearQnums(Everything& e)
