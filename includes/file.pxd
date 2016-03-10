cdef extern from "Python.h":
    ctypedef struct FILE
    FILE* stdout
    void vfprintf(FILE* f, char* s, va_list)
    void fflush(FILE*)
    FILE* fopen(char *filename, char *mode)
    void fclose(FILE*)
