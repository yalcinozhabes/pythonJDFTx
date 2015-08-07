from libcpp.vector cimport vector as stdvector #there is another global `vector` :(
from libcpp cimport bool

from includes.file cimport FILE
from includes.memory cimport shared_ptr
from includes.vector3 cimport vector3
from electronic.SpeciesInfo cimport SpeciesInfo
from electronic.IonicMinimizer cimport IonicGradient

cdef extern from "electronic/IonInfo.h" namespace "IonInfo" nogil:
    enum IonWidthMethod:
        IonWidthEcut
        IonWidthFFTbox
        IonWidthManual

cdef extern from "electronic/IonInfo.h" nogil:
    enum CoordsType:
        CoordsLattice
        CoordsCartesian

    enum coreOverlapCheck:
        additive
        vector
        none

    enum ForcesOutputCoords:
        ForcesCoordsCartesian

    cdef cppclass IonInfo:
        stdvector[ shared_ptr[SpeciesInfo] ] species
        CoordsType coordsType
        coreOverlapCheck coreOverlapCondition
        ForcesOutputCoords forcesOutputCoords
        IonWidthMethod ionWidthMethod

        IonicGradient forces
        void printPositions(FILE*)
