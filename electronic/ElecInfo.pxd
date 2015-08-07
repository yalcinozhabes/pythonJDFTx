from libcpp cimport bool
from libcpp.vector cimport vector

from core.common cimport *
from includes.vector3 cimport vector3
from includes.string cimport string


cdef extern from "electronic/ElecInfo.h" namespace "ElecInfo" nogil:
    enum FillingsUpdate:
        ConstantFillings
        FermiFillingsMix
        FermiFillingsAux
        MaximumOverlapMethod


cdef extern from "electronic/ElecInfo.h" nogil:
    enum SpinType:
        SpinNone
        SpinZ
        SpinVector
        SpinOrbit

    cdef cppclass QuantumNumber:
        vector3[double] k
        int spin
        double weight

        QuantumNumber()
        int index()

    cdef inline vector3[double] getCoord(QuantumNumber& qnum)

    cdef cppclass ElecInfo:
        int nBands, nStates
        int nDensities, spinWeight, qWeightSum
        int qStart, qStop
        bool isMine(int q)
        int whose(int q)
        int qStartOther(int iProc)
        int qStopOther(int iProc)

        SpinType spinType
        bool spinRestricted
        double nElectrons
        vector[QuantumNumber] qnums

        bool isNoncollinear()
        int spinorLength()

        FillingsUpdate fillingsUpdate

        double kT
        double mu

        int mixInterval
        bool subspaceRotation

        bool hasU

        #I am commenting out this line as I haven't figured out how to wrap std::tuple
        # vector[std::tuple[int, int, double]] customFillings
        string initialFillingsFilename

        ElecInfo()
        void mixFillings(vector[diagMatrix]& F, Energies& ener)
        void updateFillingsEnergies(vector[diagMatrix]& F, Energies&)

        # Fermi function utilities:
        double fermi(double mu, double eps)
        double fermiPrime(double mu, double eps)
        diagMatrix fermi(double mu, const diagMatrix& eps)
        diagMatrix fermiPrime(double mu, const diagMatrix& eps)

        # Propagate matrix gradient w.r.t F to gradient w.r.t. eps (in the basis where fillings are diagonal)
        matrix fermiGrad(double mu, diagMatrix& eps, matrix& gradF)

        # Compute number of electrons for a fermi distribution with specified eigenvalues
        double nElectronsFermi(double mu, vector[diagMatrix]& eps)

        # Find the chemical potential for which the fermi distribution with specified eigenvalues adds up to nElectrons
        double findMu(vector[diagMatrix]& eps, double nElectrons)

        # Find the best fit chemical potential (and optionally the density of states) given fillings and eigenvalues
        double fitMu(vector[diagMatrix]& F, vector[diagMatrix]& eps, double* dndmu)

        int findHOMO(int q)

        # //Parallel I/O utilities for diagMatrix/matrix array (one-per-kpoint, with nBands rows and columns unless overridden):
        # void read(std::vector<diagMatrix>&, const char *fname, int nRowsOverride=0) const;
        # void read(std::vector<matrix>&, const char *fname, int nRowsOverride=0, int nColsOverride=0) const;
        # void write(const std::vector<diagMatrix>&, const char *fname, int nRowsOverride=0) const;
        # void write(const std::vector<matrix>&, const char *fname, int nRowsOverride=0, int nColsOverride=0) const;

        # k-points:
        vector3[int] kfold
        void kpointsFold()
        void kpointsReduce()
