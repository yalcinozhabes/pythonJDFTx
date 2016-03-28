# This file is modified into JDFTCalculatorCPU.pyx and JDFTCalculatorGPU.pyx
# by setup.py. `{TARGET}` is replaced by `CPU` and `GPU` for these files.
# The following line left blank intentionally

# Author: Yalcin Ozhabes
# email: yalcinozhabes@gmail.com

DEF TARGET = '{TARGET}'

include "JDFTCalculatorBase.pxi"

cdef class JDFTCalculator{TARGET}:
    # cdef VerletParams v
    cdef Everything e "e"
    cdef MPIUtil* _mpiUtil "_mpiUtil"
    cdef pyComm comm "comm"
    cdef cmpi.MPI_Comm _comm "_comm"
    cdef FILE* _globalLog "_globalLog"
    cdef FILE* _nullLog "_nullLog"
    cdef int _nThreads "_nThreads"
    cdef IonicMinimizer* imin "imin"
    cdef bool didSetupRun "didSetupRun"


    #c level functions
    def __cinit__(self, *args,**kwargs):
        self.didSetupRun = False

        IF TARGET == 'GPU':
            gpuInit(stdout)

        self.e.cntrl.basisKdep = BasisKdep_BasisKpointDep

        #coords-type Cartesian
        self.e.iInfo.coordsType = CoordsType_CoordsCartesian

        #core-overlap-check vector
        self.e.iInfo.coreOverlapCondition = Condition_vector

        #coulomb-interaction Periodic
        self.e.coulombParams.geometry =  Geometry_Periodic

        #dump End None
        self.e.dump.insert(<const pair[DumpFrequency,DumpVariable]>
                      pair[DumpFrequency,DumpVariable](DumpFreq_End, DumpNone))

        #elec-cutoff 20
        self.e.cntrl.Ecut = 20.0
        self.e.cntrl.EcutRho = 0.0

        #elec-eigen-algo Davidson
        self.e.cntrl.elecEigenAlgo = ElecEigenDavidson

        #elec-ex-corr gga-PBE
        #This command does nothing by default, constructor takes care of initialization

        #electronic-minimize energyDiffThreshold  1e-08
        self.e.elecMinParams.energyDiffThreshold = 1e-8
        self.e.elecMinParams.nIterations = 400

        #exchange-regularization WignerSeitzTruncated
        self.e.coulombParams.exchangeRegularization = XReg_WignerSeitzTruncated

        #fluid None
        self.e.eVars.fluidParams.fluidType = FluidType_FluidNone

        #fluid-ex-corr lda-TF gga-PBE
        self.e.eVars.fluidParams.exCorr.kineticType = ExCorr_KineticTF
        self.e.eVars.fluidParams.exCorr.exCorrType = ExCorr_ExCorrLDA_PZ

        #fluid-gummel-loop 10 1.000000e-05
        self.e.cntrl.fluidGummel_nIterations = 10
        self.e.cntrl.fluidGummel_Atol = 1e-5

        #forces-output-coords Positions
        self.e.iInfo.forcesOutputCoords = ForcesOutput_Cartesian

        #ion-width 0
        self.e.iInfo.ionWidthMethod = IonWidthManual

        #ionic-minimize nIterations 0
        self.e.ionicMinParams.nIterations = 0

        #kpoint   0.000000000000   0.000000000000   0.000000000000  1.00000000000000
        # self.kpts = 'Gamma'

        # kpoint-folding 1 1 1
        for i in range(3):
            self.e.eInfo.kfold[i] = 1

        # latt-move-scale 1 1 1
        for i in range(3):
            self.e.cntrl.lattMoveScale[i] = 1.0

        # latt-scale 1 1 1
        for i in range(3):
            self.e.gInfo.lattScale[i] = 1.0

        # lcao-params -1 1e-06 0.001
        # nIteration is set to 3 by the constructor of ElecVars::LCAO()
        self.e.eVars.lcaoIter = 3
        self.e.eVars.lcaoTol = 1e-06
        self.e.eInfo.kT = 1e-3

        # pcm-variant GLSSA13

        #reorthogonalize-orbitals 20 1.5
        self.e.cntrl.overlapCheckInterval = 20
        self.e.cntrl.overlapConditionThreshold = 1.5

        # spintype no-spin --> initialized in __init__()

        # subspace-rotation-factor 30
        self.e.eVars.subspaceRotationFactor = 30.0

        # symmetries automatic no
        self.e.symm.mode = SymmetryMode_None
        self.e.symm.shouldMoveAtoms = False

        # self.imin = new IonicMinimizer(self.e)

    def __dealloc__(self):
        """finalizeSystem()"""
        del self._mpiUtil
        fclose(self._nullLog)
        if not (self.imin is NULL):
            del self.imin

    cdef inline void setGlobalsInline(self):
        """
        Set the global namespace for the current calculator.

        JDFTx has variables in global namespace which we can't replicate for
        each calculator. The only way to go is to call a function that modifies
        the global namespace everytime the scheduler decides to make progress
        with a calculator.
        """
        global mpiUtil, globalLog, nullLog, nProcsAvailable
        mpiUtil = self._mpiUtil
        nProcsAvailable = self._nThreads
        if mpiUtil.isHead():
            globalLog = self._globalLog
        else:
            globalLog = self._nullLog
        nullLog = self._nullLog

    #some getters and setters for easy access from python side
    property R:
        """Lattice vectors in Angstrom. Handles the conversion internally.

        Handles both the unit conversion and row major to column major
        conversion. It is a 3x3 matrix with row vectors being lattice vectors."""
        def __get__(self):
            out = np.zeros((3,3), dtype=np.double)
            for i in range(3):
                for j in range(3):
                    out[i,j] = <double>self.e.gInfo.R(j,i) * cBohr
            return out

        def __set__(self, double[:, :] value):
            for i in range(3):
                for j in range(3):
                    (&self.e.gInfo.R(j,i))[0] = <double>value[i,j] / cBohr

    property spin:
        """ In order to match the convention of ASE, skip 0 when counting
        1: Unpolarized calculation
        2: Spin-polarized calculation
        3: Non-collinear magnetism (supports spin-orbit coupling)
        4: Non-collinear without magnetization, to allow for spin-orbit"""
        def __get__(self):
            return self.e.eInfo.spinType + 1

        def __set__(self, value):
            self.e.eInfo.spinType = value - 1

    property dragWavefunctions:
        """Drag wavefunctions when ions are moved using
        atomic orbital projections (yes by default)."""
        def __get__(self):
            return self.e.cntrl.dragWavefunctions

        def __set__(self, value):
            self.e.cntrl.dragWavefunctions = <bool>value

    property kpts:
        """k-points interface with JDFTx

        JDFTx takes k-points by explicit coordinates with weights and a grid to
        fold the k-points over. ASE has a dft.kpoints module where various
        Brillouin-zone sampling schemes are implemented.

        kpts: list of tuples of kpoints and weights, all in cell coordinates.
        kpts = [(np.asarray(k), w) for k, w in zip(kpoints, weights)]
        """
        def __get__(self):
            kpts = []
            for qnum in self.e.eInfo.qnums:
                # kpt = np.ndarray(qnum.k[0], qnum.k[1], qnum.k[2]])
                kpt = np.asarray([ qnum.k[0], qnum.k[1], qnum.k[2] ])
                weight = qnum.weight
                kpts.append((kpt, weight))
            return kpts

        def __set__(self, kpts):
            clearQnums(self.e)
            cdef QuantumNumber qnum
            if kpts=='Gamma':
                addQnum(self.e, np.zeros(3), 1.0)
                return
            # else:
            cdef double[:] kpt
            cdef double weight
            for kpt, weight in kpts:
                addQnum(self.e, kpt, weight)

    property settings:
        """Dictionary that holds the initial settings.

        The settings are fixed after the first run of setup()."""
        def __set__(self, settings):
            if self.didSetupRun:
                raise AttributeError("Can't change settings after setup() runs.")
            if 'Ecut' in settings:
                self.e.cntrl.Ecut = <double>settings['Ecut']
            if 'kT' in settings:
                self.e.eInfo.kT = <double>settings['kT']
                self.e.eInfo.mixInterval = 0
                if settings['kT'] > 0:
                    self.e.eInfo.fillingsUpdate = FermiFillingsAux
                else:
                    self.e.eInfo.fillingsUpdate = ConstantFillings
            if 'nBands' in settings:
                self.e.eInfo.nBands = <int>settings['nBands']

        def __get__(self):
            settings = dict()
            settings['Ecut']   = self.e.cntrl.Ecut
            settings['kT']     = self.e.eInfo.kT
            settings['nBands'] = self.e.eInfo.nBands
            return settings


    #python level functions
    def __init__(self, comm = None, nThreads = None, **kwargs):
        """"""
        # set up the parallelization and global values
        if comm is None:
            self.comm = mpi.COMM_WORLD
        elif isinstance(comm, pyComm):
            self.comm = comm
        else:
            raise TypeError("comm has to be of type mpi4py.MPI.Comm")
        self._comm = self.comm.ob_mpi
        self._mpiUtil = new MPIUtil(self._comm)
        if nThreads is None:
            self._nThreads = utils.nThreads(self.comm)
        elif isinstance(nThreads, int):
            self._nThreads = nThreads
        else:
            raise TypeError("nThreads has to be of type int")

        # global variables (global in JDFTx not in pythonJDFTx)
        self._nullLog = fopen("/dev/null", "w")
        if 'log' in kwargs:
            log = kwargs['log']
        if log is True:
            self._globalLog = stdout
        elif isinstance(log, str):
            self._globalLog = fopen(strToCharStar(log), "w")
        elif (log is False) or (log is None):
            self._globalLog = self._nullLog
        else:
            raise TypeError("log must be one of True, False, None or a string")

        self.spin = 1

    def setGlobals(self):
        self.setGlobalsInline()

    def add_ion(self, atom):
        """
        """
        cdef char* symbol = strToCharStar(atom.symbol)

        cdef string id
        cdef vector3[double] pos
        id.assign(symbol)
        cdef shared_ptr[SpeciesInfo] sp
        sp = findSpecies(id, self.e)

        invR = np.linalg.inv(self.R)
        positionInLatticeCoordinates = atom.position.dot(invR)
        for i in range(3):
            pos[i] = <double>positionInLatticeCoordinates[i]

        if sp != 0:
            pass
        elif 'pseudopotential' in atom.data:
            if atom.data[pseudopotential].lower() in ['uspp', 'fhi', 'upf']:
                pspFile = _makePspPath(atom.symbol, atom.data[pseudopotential].lower())
            elif os.path.exists(atom.data[pseudopotential]):
                pspFile = atom.data[pseudopotential]
            else:
                raise ValueError("Can't find file " + atom.data[pseudopotential] +
                                 " or unknown format.")

            if pspFile.lower().endswith("uspp"):
                sp = newSpecies(id, strToCharStar(pspFile), PspFormat_Uspp)
            elif pspFile.lower().endswith("fhi"):
                sp = newSpecies(id, strToCharStar(pspFile), PspFormat_Fhi)
            elif pspFile.lower().endswith("upf"):
                sp = newSpecies(id, strToCharStar(pspFile), PspFormat_UPF)
            else:
                raise ValueError("pseudopotential format is not supported\n" +
                                  pspFile)
        else:
            pspFile = _makePspPath(atom.symbol)
            sp = newSpecies(id, strToCharStar(pspFile), PspFormat_Uspp)
            self.e.iInfo.species.push_back(sp)
        deref(sp).atpos.push_back(pos)

        cdef Species_Constraint constraint
        constraint.moveScale = 0.0
        constraint.type = Species_Constraint_None
        deref(sp).constraints.push_back(constraint)

    def setup(self):
        """
        Runs Everything.setup()
        """
        if self.didSetupRun:
            raise RuntimeError("setup() has already run once")
        self.setGlobalsInline()
        self.e.setup()
        self.didSetupRun = True

    def getIonicPositions(self):
        """
        Getter for ionic positions.

        Returns np.ndarray of ionic positions in cartesian coordinates,
        in JDFTx order.
        """
        cdef extern vector3[double] operator*(matrix3[double], vector3[double]&)
        cdef shared_ptr[SpeciesInfo] sp
        atpos = []
        for i in range(self.e.iInfo.species.size()):
            sp = self.e.iInfo.species[i]
            for j in range(deref(sp).atpos.size()):
                for k in range(3):
                    atpos.append((self.e.gInfo.R * deref(sp).atpos[j])[k])
        atpos = np.asarray(atpos, dtype = np.double)
        atpos.resize((len(atpos)/3,3))
        return atpos

    def updateIonicPositions(self, double[:,:] dpos):
        self.setGlobalsInline()
        if self.imin==NULL:
            self.imin = new IonicMinimizer(self.e)
        cdef IonicGradient d
        cdef int row = 0
        d.init(self.e.iInfo)
        for i in range(d.size()):
            for j in range(d[i].size()):
                for k in range(3):
                    d[i][j][k] = dpos[row][k]
                row += 1
        self.imin.step(d, 1.0)

    def runElecMin(self):
        self.setGlobalsInline()
        if self.imin is NULL:
            self.imin = new IonicMinimizer(self.e)
        # with nogil:
        resumeOperatorThreading()
        self.imin.minimize(self.e.ionicMinParams)

    def readTotalEnergy(self):
        cdef extern double double(EnergyComponents)
        return double(self.e.ener.E)

    def readForces(self):
        cdef extern vector3[double] operator*(matrix3[double], vector3[double]&)
        forces = []
        # cdef size_t i, j, k
        for i in range(int(self.e.iInfo.forces.size())):
            for j in range(int(self.e.iInfo.forces[i].size())):
                for k in range(3):
                    forces.append((self.e.gInfo.invRT * self.e.iInfo.forces[i][j])[k])
        return forces
