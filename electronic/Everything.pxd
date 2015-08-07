from libcpp.vector cimport vector
from libcpp cimport bool

from core.Coulomb cimport CoulombParams
from core.GridInfo cimport GridInfo
from core.MinimizeParams cimport MinimizeParams
from core.EnergyComponents cimport Energies
from electronic.Basis cimport Basis
from electronic.Control cimport Control
from electronic.Dump cimport Dump
from electronic.IonInfo cimport IonInfo
from electronic.ElecVars cimport ElecVars
from electronic.ElecInfo cimport ElecInfo
from electronic.Symmetries cimport Symmetries
# from electronic.ExCorr cimport ExCorr
# from electronic.VanDerWaals cimport VanDerWaals
from includes.memory cimport shared_ptr
from includes.vector3 cimport vector3

cdef extern from "electronic/Everything.h" nogil:
    cdef cppclass Everything:
        Everything()

        Control cntrl
        Dump dump
        GridInfo gInfo
        shared_ptr[GridInfo] gInfoWfns
        vector[Basis] basis
        IonInfo iInfo
        Symmetries symm
        # Symmetries symmUnperturbed
        # ExCorr exCorr
        # vector[shared_ptr[ExCorr] ] exCorrDiff
        # shared_ptr[ExactExchange] exx
        ElecInfo eInfo
        ElecVars eVars
        Energies ener

        MinimizeParams elecMinParams
        MinimizeParams ionicMinParams
        # MinimizeParams fluidMinParams
        # MinimizeParams latticeMinParams
        # MinimizeParams inverseKSminParams
        # VerletParams verletParams
        # SCFparams scfParams
        #
        CoulombParams coulombParams
        # shared_ptr[Coulomb] coulomb
        #
        # shared_ptr[VanDerWaals] vanDerWaals

        void setup() except +
        void updateSupercell(bool force=false)
