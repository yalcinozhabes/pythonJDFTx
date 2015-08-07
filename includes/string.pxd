cdef extern from "<string>" namespace "std" nogil:
    cdef cppclass basic_string[_Char,_Traits]:
        basic_string() except +
        basic_string& assign(char* ) except +
        basic_string& operator+(char*) except +

cdef extern from "core/string.h" nogil:
    cdef struct ichar_traits
    ctypedef basic_string[char, ichar_traits] string
