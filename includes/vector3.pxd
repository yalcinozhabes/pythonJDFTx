cdef extern from "core/vector3.h":
    cdef cppclass vector3[scalar]:
#        scalar v[3]
#
        vector3()
        scalar& operator[](int)
        scalar& operator[](size_t)
