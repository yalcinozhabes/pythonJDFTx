from libcpp cimport bool
from includes.file cimport FILE

cdef extern from "core/GpuUtil.h" nogil:
    bool gpuInit(FILE *)
