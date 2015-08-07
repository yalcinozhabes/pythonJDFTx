from libcpp.map cimport map
from includes.string cimport string

cdef extern from "core/EnergyComponents.h" nogil:
    cdef cppclass EnergyComponents(map[string, double]):
        pass

cdef extern from "electronic/Energies.h" nogil:
    cdef struct Energies:
        EnergyComponents E
