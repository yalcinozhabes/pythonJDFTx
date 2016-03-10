# Author: Yalcin Ozhabes
# email: yalcinozhabes@gmail.com
include "electronic/QuantumNumber.pyx"

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
from electronic.SpeciesInfo cimport Uspp as PseudopotentialFormat_Uspp
from electronic.SpeciesInfo cimport Constraint as Species_Constraint
from electronic.SpeciesInfo cimport None as Species_Constraint_None
from electronic.Dump cimport DumpFrequency, DumpFreq_End, DumpVariable, DumpNone
from electronic.QuantumNumber cimport QuantumNumber
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

cimport numpy as np
import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree

#globals
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
    deref(specie).pspFormat = PseudopotentialFormat_Uspp
    return specie

from cpython cimport PY_MAJOR_VERSION
cdef char* strToCharStar(object text):
    if isinstance(text, unicode): # most common case first
        utf8_data = text.encode('UTF-8')
    elif (PY_MAJOR_VERSION < 3) and isinstance(text, str):
        text.decode('ASCII') # trial decoding, or however you want to check for plain ASCII data
        utf8_data = text
    else:
        raise ValueError("requires text input, got %s" % type(text))
    return utf8_data

cdef class JDFTCalculator:
    # cdef VerletParams v
    cdef Everything e
    cdef MPIUtil* _mpiUtil
    cdef FILE* _globalLog
    cdef FILE* _nullLog
    cdef IonicMinimizer* imin

    #c level functions
    def __cinit__(self, *args,**kwargs):
        self._nullLog = fopen("/dev/null", "w")
        self._globalLog = stdout

        self._mpiUtil = new MPIUtil(mpi.MPI_COMM_WORLD) #initSystemCmdLine

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
        global mpiUtil, globalLog, nullLog
        mpiUtil = self._mpiUtil
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

    #python level functions
    def __init__(self, *args, **kwargs):
        self.spin = 1

    def setGlobals(self):
        self.setGlobalsInline()

    def disableLog(self):
        self._globalLog = self._nullLog

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
        positionInLatticeCoordinates = invR.dot(atom.position)
        for i in range(3):
            pos[i] = <double>positionInLatticeCoordinates[i]

        if sp != 0:
            deref(sp).atpos.push_back(pos)
        else:
            pspFile = _makePspPath(atom.symbol)
            sp = createNewSpecie(id, strToCharStar(pspFile))
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
        self.setGlobalsInline()
        print("globals are set")
        self.e.setup()

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
