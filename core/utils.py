# email: yalcinozhabes@gmail.com
# Author: Yalcin Ozhabes

from __future__ import division

import socket, multiprocessing
from mpi4py import MPI as mpi


def nThreads(comm):
    """
    Return how many cpu cores are available
    for this particular process

    Divide number of total cores to number of processes on this node.
    """
    hostname = socket.gethostname()
    hostnames = ['' for _ in range(comm.size)]

    for r in range(comm.size):
        if comm.rank == r:
            hostnames[r] = hostname
        hostnames[r] = comm.bcast(hostname, root=r)
    N = multiprocessing.cpu_count()
    nSiblings = 0 # total number of processes on this node
    iSibling = 0 # index of the process among the ones on this node
    for r, host in enumerate(hostnames):
        if host==hostname:
            if comm.rank == r:
                iSibling = nSiblings
            nSiblings += 1
    return N*(iSibling+1)//nSiblings - N*iSibling//nSiblings
