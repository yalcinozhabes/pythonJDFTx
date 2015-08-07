cimport ManagedMemory
from cython.operator cimport dereference as deref
from libcpp cimport bool
from libc.stdlib cimport malloc, free
from log cimport *

cdef class ManagedMemory:
    thisptr = ManagedMemory()
    def __cinit__(self,*args,**kwargs):
        cdef int argc = len(args)
        cdef char** argv = <char**>malloc(argc * sizeof(char*))

        for i in range(argc):
            argv[i] = args[i]

        initSystemCmdline(
                argc, argv,
                "Performs Joint Density Functional Theory calculations.",
                self.inputFilename, self.dryRun, self.printDefaults
                )
        parse(readInputFile(self.inputFilename), self.e, self.printDefaults);
        self.e.setup()

        if self.dryRun:
            logPrintf("Dry run successful: commands are valid and initialization succeeded.\n");
        elif self.e.cntrl.dumpOnly:
            logPrintf("\n----------- Energy evaluation at fixed state -------------\n")
            logFlush()
        elif self.e.cntrl.fixed_H:
            pass
        elif self.e.vibrations:
            pass
        elif self.e.latticeMinParams.nIterations:
            pass
        elif self.e.verletParams.tMax:
            pass
        else:
            imin = new IonicMinimizer(self.e)
            imin.minimize(self.e.ionicMinParams)
        self.e.dump(DumpFreq_End, 0)
        finalizeSystem()
        return

    def __init__(self,*args,**kwargs):
        pass
