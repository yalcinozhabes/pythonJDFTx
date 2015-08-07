from libcpp.vector cimport vector

from core.MinimizeParams cimport MinimizeParams
from electronic.Everything cimport Everything
from electronic.IonInfo cimport IonInfo
from includes.vector3 cimport vector3

cdef extern from "electronic/IonicMinimizer.h" nogil:
    cdef cppclass IonicGradient(vector[vector[vector3[double]]]):
        void init(IonInfo&)

    cdef cppclass IonicMinimizer:
        IonicMinimizer(Everything& e)
        double minimize(MinimizeParams& params)
        void step(IonicGradient& dir, double alpha);
