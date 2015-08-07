from libcpp.vector cimport vector
from libcpp cimport bool

from includes.memory cimport shared_ptr
from includes.string cimport string
from electronic.ExCorr cimport ExCorr
from fluid.FluidComponent cimport FluidComponent

cdef extern from "fluid/FluidSolverParams.h" nogil:
    enum FluidType:
        FluidNone
        FluidLinearPCM
        FluidNonlinearPCM
        FluidSaLSA
        FluidClassicalDFT

    enum FMixFunctional:
        FMixNone
        LJPotential
        GaussianKernel

    cdef struct FluidSolverParams:
        FluidType fluidType
        # PCMVariant pcmVariant

        double T
        double P
        double epsBulkOverride, epsInfOverride
        bool verboseLog

        vector[ shared_ptr[FluidComponent] ]& components
        vector[ shared_ptr[FluidComponent] ]& solvents
        vector[ shared_ptr[FluidComponent] ]& cations
        vector[ shared_ptr[FluidComponent] ]& anions

        void addComponent(const shared_ptr[FluidComponent]& component)

        #Fit parameters:
        double nc
        double sigma
        double cavityTension
        double vdwScale

        #For CANDLE alone:
        double Ztot
        double eta_wDiel
        double sqrtC6eff
        double pCavity

        #For SCCS alone:
        double rhoMin, rhoMax
        double rhoDelta
        double cavityPressure

        #For SaLSA alone:
        int lMax

        #Debug parameters for Nonlinear PCM's:
        bool linearDielectric
        bool linearScreening
        bool nonlinearSCF
        #PulayParams scfParams

        #For Explicit Fluid JDFT alone:
        ExCorr exCorr
        #vector[FmixParams] FmixList

        string initWarnings

        FluidSolverParams()
        void setPCMparams()
        void setCDFTparams()
        bool needsVDW()
        bool ionicScreening()
