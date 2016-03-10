from libcpp cimport bool
cimport mpi4py.libmpi as mpi

cdef extern from "core/MPIUtil.h" nogil:
    cdef cppclass MPIUtil:
        MPIUtil(mpi.MPI_Comm)
        bool isHead()
        int nProcesses() const
