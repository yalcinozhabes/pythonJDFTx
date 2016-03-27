# email: yalcinozhabes@gmail.com
# Author: Yalcin Ozhabes
include "electronic/QuantumNumber.pyx"

from cython.operator cimport dereference as deref
from libcpp cimport bool
from libcpp.pair cimport pair
from mpi4py.MPI cimport Comm as pyComm
cimport mpi4py.libmpi as cmpi
cdef extern from "includes/mpi-compat.h": pass

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
from core.MPIUtil cimport MPIUtil
from core.Thread cimport *
from core.Util cimport *
IF TARGET=='GPU':
    from core.GpuUtil cimport gpuInit
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
import multiprocessing, socket
from mpi4py import MPI as mpi

from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree

from core import utils

#globals
cdef:
    double cBohr "Bohr" = Bohr
    double cHartree "Hartree" = Hartree

def _makePspPath(symbol):
    import os
    pspFolder = os.getenv("PSEUDOPOT_HOME")
    return os.path.join(pspFolder, symbol.lower() + ".uspp")

cdef shared_ptr[SpeciesInfo] createNewSpecies "createNewSpecies" (string id, char* pspFile):
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
