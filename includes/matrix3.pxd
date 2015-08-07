from includes.vector3 cimport vector3

cdef extern from "core/matrix3.h":
    cdef cppclass matrix3[scalar]:
        scalar m[3][3]
        matrix3()
        matrix3(scalar, scalar, scalar)
        scalar& operator()(int, int)
        
