from libcpp cimport bool

cdef extern from "<WBCKits.hpp>" namespace "WBCKits":

    cdef enum AlgorithmOfHQP:
        HQP_Origin
        HQP_Nullspace
    
    cdef enum DecomposeMethod:
        Decompose_QR = 0
        Decompose_SVD
        Decompose_COD

    cdef cppclass Options:
        Options()
        Options(const Options&)
        # Options& operator=( const Options&)  # equality operator cannot be overloaded in Python
        bool SetDefault()
        void Print() const

        # public values
        int nwsr_max_
        double cpu_time_max_

        bool enable_regulation_
        double weight_regulation_

        AlgorithmOfHQP hqp_algorithm_
        DecomposeMethod decompose_method_
        double decompose_threshold_

        bool get_slack_variable_

    cdef cppclass WBCBase:
        # WBCBase()
        # bool SolveWBC()
        # void AddTask(const MatrixXd &A, const VectorXd &b)
        # void AddTask(const MatrixXd &A, const VectorXd &b, const VectorXd &w)
        # void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC, bool hard_constraint = true)
        # void Reset()
        pass

    cdef cppclass WQP:
        WQP()
        void AddTask(const MatrixXd &A, const VectorXd &b)
        void AddTask(const MatrixXd &A, const VectorXd &b, const VectorXd &w)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC, bool hard_constraint)
        void Reset()
        bool SolveWBC()
        void GetResult(VectorXd & optimal_result)
        bool SetOptions( Options& options)
        double GetCostOfQP() const

    cdef cppclass RHP:
        RHP()
        void AddTask(const MatrixXd &A, const VectorXd &b)
        void AddTask(const MatrixXd &A, const VectorXd &b, const VectorXd &w)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC, bool hard_constraint)
        void Reset()
        void SetPriority(const MatrixXd &phi, const VectorXi &chi)
        void SetPriority(const VectorXi &phi, const VectorXi &chi)
        bool SolveWBC()
        void GetResult(VectorXd & optimal_result)
        bool SetOptions( const Options& options)
        double GetCostOfQP() const

    cdef cppclass HQP:
        HQP()
        void AddTask(const MatrixXd &A, const VectorXd &b)
        void AddTask(const MatrixXd &A, const VectorXd &b, const VectorXd &w)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC, bool hard_constraint)
        void Reset()
        void SetPriority(const VectorXi &phi, const VectorXi &chi)
        bool SolveWBC()
        void GetResult(VectorXd & optimal_result)
        bool SetOptions( const Options& options)
        double GetCostOfQP() const

    cdef cppclass GHC:
        GHC()
        void AddTask(const MatrixXd &A, const VectorXd &b)
        void AddTask(const MatrixXd &A, const VectorXd &b, const VectorXd &w)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC)
        void AddConstraint(const MatrixXd &C, const VectorXd &lbC, const VectorXd &ubC, bool hard_constraint)
        void Reset()
        void SetPriority(const MatrixXd &phi, const VectorXi &chi)
        void SetPriority(const VectorXi &phi, const VectorXi &chi)
        bool SolveWBC()
        void GetResult(VectorXd & optimal_result)
        bool SetOptions( const Options& options)
        double GetCostOfQP() const

cdef extern from "<eigen3/Eigen/Dense>" namespace "Eigen":
    cdef cppclass MatrixXd:
        MatrixXd()
        MatrixXd (int rows, int cols)
        int rows()
        int cols()
        void resize (int,int)
        double& coeff "operator()"(int,int)
        double* data()
        void setZero()

    cdef cppclass VectorXd:
        VectorXd()
        VectorXd (int dim)
        int rows()
        int cols()
        void resize (int)
        double& operator[](int)
        double* data()
    
    cdef cppclass VectorXi:
        VectorXi()
        VectorXi (int dim)
        int rows()
        int cols()
        void resize (int)
        double& operator[](int)
        double* data()
