from libcpp.vector cimport vector
from libcpp cimport bool

from includes.vector3 cimport vector3
from includes.matrix3 cimport matrix3
from includes.memory cimport shared_ptr
from core.common cimport *

#include <core/ScalarField.h]
#include <core/GpuUtil.h>
#include <map>
#include <fftw3.h>
#include <stdint.h>
#include <cstdio>
#include <mutex>

cdef extern from "core/GridInfo.h" nogil:
    cdef cppclass GridInfo:
        vector3[double] lattScale
        vector3[int] S
        matrix3[double] R
        double Gmax
        double GmaxRho
        double detR
        matrix3[double] RT, RTR, invR, invRT, invRTR
        matrix3[double] G, GT, GGT, invGGT
        double dV
        vector3[double] h[3]
        int nr
        int nG
        double dGradial
        double GmaxSphere
        double GmaxGrid
        double maxAllowedStrain

        GridInfo()

        void update()
        void setLatticeVectors()
        void initialize(bool skipHeader=false, const vector[ matrix3[int] ] sym = vector[ matrix3[int] ](1, matrix3[int](1,1,1)))

        int irStart, irStop
        int iGstart, iGstop

        inline vector3[int] wrapGcoords(const vector3[int] iG) const
        inline int fullRindex(const vector3[int] iR) const
        inline int fullGindex(const vector3[int] iG) const
        inline int halfGindex(const vector3[int] iG) const
