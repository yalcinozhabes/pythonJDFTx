from libcpp.vector cimport vector

from includes.vector3 cimport vector3
from core.GridInfo cimport GridInfo
from electronic.IonInfo cimport IonInfo


cdef extern from "electronic/Basis.h" nogil:
    cdef cppclass Basis:
        const GridInfo* gInfo
        const IonInfo* iInfo

        size_t nbasis
        vector3[int] *iGarr
        int *index
        int *indexPref
        vector3[int] *iGarrPref

        vector[int] head

        Basis()
        Basis(const Basis&)

        void setup(const GridInfo& gInfo, const IonInfo& iInfo, double Ecut, const vector3[double] k)
        void setup(const GridInfo& gInfo, const IonInfo& iInfo, const vector[int]& indexVec)
