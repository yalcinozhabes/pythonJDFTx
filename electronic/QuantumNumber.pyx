cimport electronic.QuantumNumber
from electronic.Everything cimport Everything
cimport numpy as cnp

cdef void addQnum(Everything& e, double[:] kpt, double weight):
    """add kpoint with weight.

    kpt must be an ndarray of size 3 representing a vector in cell coordinates.
    """
    cdef QuantumNumber qnum
    cdef int i
    for i in range(3):
        (&(qnum.k[i]))[0] = kpt[i]
        qnum.weight = weight
    e.eInfo.qnums.push_back(qnum)
    e.eInfo.nStates = e.eInfo.qnums.size()

cdef void clearQnums(Everything& e):
    e.eInfo.qnums.clear()
    e.eInfo.nStates = 0
