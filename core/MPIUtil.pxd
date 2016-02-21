cimport mpi4py.libmpi as mpi

cdef extern from "core/MPIUtil.h" nogil:
    cdef cppclass MPIUtil:
        MPIUtil(mpi.MPI_Comm)
        int nProcesses() const
