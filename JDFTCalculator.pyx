# Author: Yalcin Ozhabes
# email: yalcinozhabes@gmail.com

from cython.operator cimport dereference as deref
from libcpp cimport bool
from libcpp.pair cimport pair
cimport mpi4py.libmpi as mpi

from includes.string cimport string
from includes.file cimport *
from includes.memory cimport shared_ptr
from includes.vector3 cimport vector3
from includes.matrix3 cimport matrix3
from electronic.Everything cimport Everything
from electronic.Basis cimport Basis
from electronic.SpeciesInfo cimport SpeciesInfo, findSpecies
from electronic.SpeciesInfo cimport Uspp as PseudoPotentialFormat_Uspp
from electronic.SpeciesInfo cimport Constraint as Species_Constraint
from electronic.SpeciesInfo cimport None as Species_Constraint_None
from electronic.Dump cimport DumpFrequency, DumpFreq_End, DumpVariable, DumpNone
from electronic.ElecInfo cimport QuantumNumber
# from electronic.ExCorr cimport ExCorr
from core.MPIUtil cimport MPIUtil
from core.Thread cimport *
from core.Util cimport *
from core.GridInfo cimport GridInfo
from core.Coulomb cimport CoulombParams
from core.Coulomb cimport Periodic as Geometry_Periodic
from core.Coulomb cimport WignerSeitzTruncated as XReg_WignerSeitzTruncated
from core.EnergyComponents cimport EnergyComponents
from electronic.Control cimport Control, ElecEigenDavidson
from electronic.Control cimport BasisKpointDep as BasisKdep_BasisKpointDep
from electronic.IonInfo cimport IonInfo
from electronic.IonInfo cimport CoordsCartesian as CoordsType_CoordsCartesian
from electronic.IonInfo cimport ForcesCoordsCartesian as ForcesOutput_Cartesian
from electronic.IonInfo cimport vector as Condition_vector
from electronic.IonInfo cimport IonWidthManual
from electronic.ExCorr cimport KineticTF as ExCorr_KineticTF
from electronic.ExCorr cimport ExCorrLDA_PZ as ExCorr_ExCorrLDA_PZ
from electronic.IonicMinimizer cimport IonicMinimizer, IonicGradient
from electronic.Symmetries cimport SymmetriesNone as SymmetryMode_None
from fluid.FluidParams cimport FluidSolverParams
from fluid.FluidParams cimport FluidNone as FluidType_FluidNone

import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree
cdef:
    double cBohr = Bohr
    double cHartree = Hartree
    int i, j, k


def _makePspPath(symbol):
    import os
    pspFolder = os.getenv("PSEUDOPOT_HOME")
    return os.path.join(pspFolder, symbol.lower() + ".uspp")

cdef shared_ptr[SpeciesInfo] createNewSpecie(string id, char* pspFile):
    # cdef SpeciesInfo* speciePtr = new SpeciesInfo()
    cdef shared_ptr[SpeciesInfo] specie = shared_ptr[SpeciesInfo](new SpeciesInfo())
    deref(specie).potfilename.assign(pspFile)
    deref(specie).fromWildcard = False
    deref(specie).name = id
    deref(specie).pspFormat = PseudoPotentialFormat_Uspp
    return specie

cdef class JDFTCalculator:
    # cdef VerletParams v
    cdef Everything e
    cdef IonicMinimizer* imin

    #c level functions
    def __cinit__(self, *args,**kwargs):
        global mpiUtil, globalLog, nullLog
        resumeOperatorThreading()

        nullLog = fopen("/dev/null","w")
        # globalLog = nullLog
        mpiUtil = new MPIUtil(mpi.MPI_COMM_WORLD) #initSystemCmdLine

        #Default commands that sets some variables up at the beginning:
        #basis kpoint-dependent
        self.e.cntrl.basisKdep = BasisKdep_BasisKpointDep

        #coords-type Cartesian
        self.e.iInfo.coordsType = CoordsType_CoordsCartesian

        #core-overlap-check vector
        self.e.iInfo.coreOverlapCondition = Condition_vector

        #coulomb-interaction Periodic
        self.e.coulombParams.geometry =  Geometry_Periodic

        #dump End None
        self.e.dump.insert(pair[DumpFrequency,DumpVariable](DumpFreq_End, DumpNone))

        #elec-cutoff 20
        self.e.cntrl.Ecut = 20.0
        self.e.cntrl.EcutRho = 0.0

        #elec-eigen-algo Davidson
        self.e.cntrl.elecEigenAlgo = ElecEigenDavidson

        #elec-ex-corr gga-PBE
        #This command does nothing by default, constructor takes care of initialization

        #electronic-minimize energyDiffThreshold  1e-08
        self.e.elecMinParams.energyDiffThreshold = 1e-8

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
        cdef QuantumNumber qnum
        for i in range(3):
            (&(qnum.k[i]))[0] = 0.0
            qnum.weight = 1.0
        self.e.eInfo.qnums.push_back(qnum)
        self.e.eInfo.nStates = self.e.eInfo.qnums.size()

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
        global mpiUtil
        del mpiUtil
        if not (self.imin is NULL):
            del self.imin

    #some getters and setters for easy access from python side
    property R:
        """Lattice vectors in Angstrom. Handles the conversion internally."""
        def __get__(self):
            cdef double unitConvertedEntry
            out = np.zeros((3,3), dtype=np.double)
            for i in range(3):
                for j in range(3):
                    unitConvertedEntry = self.e.gInfo.R(i,j) * cBohr
                    out[i,j] = unitConvertedEntry
            return out

        def __set__(self, value):
            cdef double unitConvertedEntry
            cdef double& tmpR
            for i in range(3):
                for j in range(3):
                    unitConvertedEntry = value[i,j] / cBohr
                    # self.e.gInfo.R(i,j) = unitConvertedEntry
                    (&self.e.gInfo.R(i,j))[0] = unitConvertedEntry

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

    #python level functions
    def __init__(self, *args, **kwargs):
        self.spin = 1

    def disableLog(self):
        global globalLog, nullLog
        globalLog = nullLog

    def add_ion(self, atom):
        symbol = str(atom.symbol)

        cdef string id
        cdef vector3[double] pos
        id.assign(<char*>symbol)
        cdef shared_ptr[SpeciesInfo] sp
        sp = findSpecies(id, self.e)

        invR = np.linalg.inv(self.R)
        positionInLatticeCoordinates = invR.dot(atom.position)
        for i in range(3):
            pos[i] = <double>positionInLatticeCoordinates[i]

        if sp != 0:
            deref(sp).atpos.push_back(pos)
        else:
            pspFile = _makePspPath(atom.symbol)
            sp = createNewSpecie(id, <char*>pspFile)
            self.e.iInfo.species.push_back(sp)
            deref(sp).atpos.push_back(pos)

        cdef Species_Constraint constraint
        constraint.moveScale = 0.0
        constraint.type = Species_Constraint_None
        deref(sp).constraints.push_back(constraint)

    def setup(self):
        """Runs Everything.setup"""
        self.e.setup()

    def readIonicPositions(self):
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

    def updateIonicPositions(self,double[:,:] dpos):
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
