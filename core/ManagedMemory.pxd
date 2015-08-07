cdef extern from "core/ColumnBundle.h" nogil:
    cdef cppclass ColumnBundle:
        complex* data()
        const complex* data() const

        size_t nData() const
        bool isOnGpu() const

        #ifdef GPU_ENABLED
        complex* dataGpu()
        const complex* dataGpu() const
        #endif

        inline complex* dataPref()
        inline const complex* dataPref() const

        void send(int dest, int tag=0) const
        void recv(int src, int tag=0)
        void bcast(int root=0)
        void allReduce(MPIUtil::ReduceOp op, bool safeMode=false, bool ignoreComplexCheck=false)

        void write(const char *fname) const
        void writea(const char *fname) const
        void write(FILE *filep) const
        void read(const char *fname)
        void read(FILE *filep)
        void read_real(const char *fname)
        void read_real(FILE *filep)
        void write_real(const char *fname) const
        void write_real(FILE *filep) const
        void dump(const char* fname, bool realPartOnly) const
        void zero()

        static void reportUsage()

        int nCols() const 
        size_t colLength() const
        explicit operator bool() const

        size_t index(int i, size_t j) const

        bool isSpinor() const
        int spinorLength() const

        const QuantumNumber *qnum
        const Basis *basis

        void init(int nc, size_t len, const Basis* b, const QuantumNumber* q, bool onGpu=false)
        void free()
        ColumnBundle(int nc=0, size_t len=0, const Basis* b=NULL, const QuantumNumber* q=NULL, bool onGpu=false)
        ColumnBundle(const ColumnBundle&)
        ColumnBundle(ColumnBundle&&)

        ColumnBundle similar(int ncOverride=-1) const

        ColumnBundle& operator=(const ColumnBundle&)
        ColumnBundle& operator=(ColumnBundle&&)

        ColumnBundle getSub(int colStart, int colStop) const
        void setSub(int colStart, const ColumnBundle&)

        complexScalarFieldTilde getColumn(int i, int s) const
        void setColumn(int i, int s, const complexScalarFieldTilde&)
        void accumColumn(int i, int s, const complexScalarFieldTilde&)

        void randomize(int colStart, int colStop)
