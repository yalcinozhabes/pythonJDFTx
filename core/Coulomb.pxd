cdef extern from "core/Coulomb.h" namespace "CoulombParams" nogil:
    enum Geometry:
        Periodic
        Slab
        Wire
        Cylindrical
        Isolated
        Spherical
    enum ExchangeRegularization:
        WignerSeitzTruncated

cdef extern from "core/Coulomb.h" nogil:
    cdef struct CoulombParams:
        Geometry geometry
        ExchangeRegularization exchangeRegularization
