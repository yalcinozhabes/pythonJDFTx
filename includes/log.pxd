#Utility functions to use C++ level print functions
#DEPRECATED

cdef extern from "stdarg.h":
    ctypedef struct va_list:
        pass
    ctypedef struct fake_type:
        pass
    void va_start(va_list, void* arg)
    void* va_arg(va_list, fake_type)
    void va_end(va_list)
    fake_type int_type "int"

cdef extern from "Python.h":
    ctypedef struct FILE
    void vfprintf(FILE* f, char* s, va_list)
    void fflush(FILE*)

cdef extern from "core/Util.h":
    FILE* globalLog

cdef inline void logPrintf(char* s, ...):
    cdef va_list args;
    va_start(args, s);
    vfprintf(globalLog, s, args)
    va_end(args)
    return
cdef inline void logFlush():
    fflush(globalLog)
