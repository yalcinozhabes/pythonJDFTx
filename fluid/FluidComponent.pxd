from includes.memory cimport shared_ptr

cdef extern from "fluid/FluidComponent.h" namespace "FluidComponent" nogil:
    enum Name:
        #Neutral solvent molecules:
        H2O
        CHCl3
        CCl4
        CH3CN
        DMC
        EC
        PC
        DMF
        THF
        EthylEther
        Chlorobenzene
        Isobutanol
        CarbonDisulfide
        DMSO
        CH2Cl2
        Ethanol
        Methanol
        Octanol
        Glyme
        EthyleneGlycol
        CustomSolvent
        #Cations:
        Sodium
        HydratedSodium
        Potassium
        HydratedPotassium
        Hydronium
        HydratedHydronium
        CustomCation
        #Anions:
        Chloride
        Fluoride
        Perchlorate
        Hydroxide
        HydratedHydroxide
        CustomAnion

    enum Type:
        Solvent
        Cation
        Anion

    enum Functional:
        ScalarEOS
        FittedCorrelations
        BondedVoids
        MeanFieldLJ
        FunctionalNone

    enum Representation:
        Pomega
        PsiAlpha
        MuEps

    enum TranslationMode:
        ConstantSpline
        LinearSpline
        Fourier


cdef extern from "fluid/FluidComponent.h" nogil:
    cdef cppclass FluidComponent:
        Name name
        Type type
        Type getType(Name name)

        const Functional functional
        double epsLJ

        Representation representation

        TranslationMode translationMode

        #Bulk solvent properties (used by various PCM's):
        double epsBulk
        double Nbulk
        double pMol
        double epsInf
        double Pvap
        double sigmaBulk
        double Rvdw
        double Res

        double Nnorm

        #Molecule geometry and site properties:
        # Molecule molecule

        double pureNbulk(double T)

        FluidComponent(Name name, double T, Functional functional)

        #Extra properties when participating in a classical density functional FluidMixture:
        # shared_ptr[class SO3quad] quad
        # shared_ptr[class TranslationOperator] trans
        # shared_ptr[class IdealGas] idealGas
        # shared_ptr[class Fex] fex
        # shared_ptr[class ScalarEOS] eos
        # unsigned offsetIndep
        # unsigned offsetDensity
        # void addToFluidMixture(class FluidMixture* fluidMixture)
