from libcpp.vector cimport vector
from libcpp cimport bool

from includes.memory cimport shared_ptr
from includes cimport file
from includes.matrix3 cimport matrix3
from core.common cimport *
from electronic.Everything cimport Everything
from electronic.Basis cimport Basis

cdef extern from "electronic/Symmetries.h" nogil:
    enum SymmetryMode:
        SymmetriesNone
        SymmetriesAutomatic
        SymmetriesManual

    cdef cppclass Symmetries:
        SymmetryMode mode
        bool shouldMoveAtoms

        Symmetries()
        void setup(const Everything& everything)
        void setupMesh()


        vector[QuantumNumber] reduceKmesh(const vector[QuantumNumber]& qnums) const

        # void symmetrize(ScalarField&) const
        # void symmetrize(IonicGradient&) const
        # void symmetrizeSpherical(matrix&, const SpeciesInfo* specie) const
        # const vector[ matrix3[int] ]& getMatrices() const
        # const vector[ matrix3[int] ]& getMeshMatrices() const
        # const vector[matrix]& getSphericalMatrices(int l, bool relativistic) const
        # const vector[int]& getKpointInvertList() const
        # const vector[vector[vector[int] ] ]& getAtomMap() const
        # void printKmap(FILE* fp) const
