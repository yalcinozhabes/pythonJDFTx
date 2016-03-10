from libcpp.vector cimport vector
from libcpp cimport bool

from includes.memory cimport shared_ptr
from includes.string cimport string
from includes.vector3 cimport vector3
from core.common cimport *
from electronic.Everything cimport Everything

cdef extern from "electronic/SpeciesInfo.h" namespace "SpeciesInfo::Constraint" nogil:
    enum ConstraintType:
        None
        Linear
        Planar

cdef extern from "electronic/SpeciesInfo.h" namespace "SpeciesInfo" nogil:
    enum PseudopotentialFormat:
        Fhi
        Uspp
        UPF
    struct Constraint:
        double moveScale
        ConstraintType type


cdef extern from "electronic/SpeciesInfo.h" nogil:
    cdef cppclass SpeciesInfo:
        string name
        string potfilename, pulayfilename
        bool fromWildcard
        PseudopotentialFormat pspFormat

        vector[vector3[double]] atpos
        vector[Constraint] constraints
        void setup(const Everything&);

cdef extern from "commands/command.h" nogil:
    cdef shared_ptr[SpeciesInfo] findSpecies(string id, Everything& e)
