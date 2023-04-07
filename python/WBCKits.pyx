# distutils: language = c++

import numpy as np
cimport numpy as np
from libcpp cimport bool

from cython.operator cimport dereference as deref

cimport cWBCKits

##############################
# Conversion Numpy <-> Eigen
##############################

# VectorXd
cdef cWBCKits.VectorXd NumpyToVectorXd (np.ndarray[np.double_t, ndim=1, mode="c"] x):
    cdef cWBCKits.VectorXd cx = cWBCKits.VectorXd(x.shape[0])
    for i in range (x.shape[0]):
        cx[i] = x[i]

    return cx

cdef np.ndarray VectorXdToNumpy (cWBCKits.VectorXd cx):
    result = np.ndarray ((cx.rows()))
    for i in range (cx.rows()):
        result[i] = cx[i]

    return result

# MatrixXd
cdef cWBCKits.MatrixXd NumpyToMatrixXd (np.ndarray[np.double_t, ndim=2, mode="c"] M):
    cdef cWBCKits.MatrixXd cM = cWBCKits.MatrixXd(M.shape[0], M.shape[1])
    for i in range (M.shape[0]):
        for j in range (M.shape[1]):
            (&(cM.coeff(i,j)))[0] = M[i,j]

    return cM

# cdef np.ndarray MatrixXdToNumpy (cWBCKits.MatrixXd cM):
#     result = np.ndarray ([cM.rows(), cM.cols()])
#     for i in range (cM.rows()):
#         for j in range (cM.cols()):
#             result[i,j] = cM.coeff(i,j)

#     return result

# VectorXi

##############################
#  enumerate types wrapper
##############################
# cdef class PyAlgorithmOfHQP:
HQP_Origin      = cWBCKits.HQP_Origin
HQP_Nullspace   = cWBCKits.HQP_Nullspace

# cdef class pyDecomposeMethod:
Decompose_QR    = cWBCKits.Decompose_QR
Decompose_SVD   = cWBCKits.Decompose_SVD
Decompose_COD   = cWBCKits.Decompose_COD

##############################
#  WBCKits wrapper
##############################
cdef class PyOptions:
    cdef cWBCKits.Options *thisptr
    def __cinit__(self):
        self.thisptr = new cWBCKits.Options()
    def __dealloc__(self):
        del self.thisptr
    
    def SetDefault(self):
        return self.thisptr.SetDefault()
    
    def Print(self):
        return self.thisptr.Print()

    # public values
    property nwsr_max_:
        def __get__(self): return self.thisptr.nwsr_max_
        def __set__(self, nwsr_max_): self.thisptr.nwsr_max_ = nwsr_max_
    
    property cpu_time_max_:
        def __get__(self): return self.thisptr.cpu_time_max_
        def __set__(self, cpu_time_max_): self.thisptr.cpu_time_max_ = cpu_time_max_

    property enable_regulation_:
        def __get__(self): return self.thisptr.enable_regulation_
        def __set__(self, enable_regulation_): self.thisptr.enable_regulation_ = enable_regulation_

    property weight_regulation_:
        def __get__(self): return self.thisptr.weight_regulation_
        def __set__(self, weight_regulation_): self.thisptr.weight_regulation_ = weight_regulation_

    property hqp_algorithm_:
        def __get__(self): return self.thisptr.hqp_algorithm_
        def __set__(self, hqp_algorithm_): self.thisptr.hqp_algorithm_ = hqp_algorithm_

    property decompose_method_:
        def __get__(self): return self.thisptr.decompose_method_
        def __set__(self, decompose_method_): self.thisptr.decompose_method_ = decompose_method_

    property decompose_threshold_:
        def __get__(self): return self.thisptr.decompose_threshold_
        def __set__(self, decompose_threshold_): self.thisptr.decompose_threshold_ = decompose_threshold_
    
    property get_slack_variable_:
        def __get__(self): return self.thisptr.get_slack_variable_
        def __set__(self, get_slack_variable_): self.thisptr.get_slack_variable_ = get_slack_variable_


cdef class WBCBase:
    pass

cdef class WQP:
    cdef cWBCKits.WQP *thisptr
    # construct function in python
    def __cinit__(self):
        self.thisptr = new cWBCKits.WQP()

    # destruct function in python
    def __dealloc__(self):
        del self.thisptr
    # declare class methods    
    def AddTask(self,
                np.ndarray[np.double_t, ndim=2, mode="c"] A,
                np.ndarray[np.double_t, ndim=1, mode="c"] b,
                np.ndarray[np.double_t, ndim=1, mode="c"] w=None):
        if w is None:
            self.thisptr.AddTask( NumpyToMatrixXd(A),
                                  NumpyToVectorXd(b) )
        else:
            self.thisptr.AddTask( NumpyToMatrixXd(A),
                                  NumpyToVectorXd(b),
                                  NumpyToVectorXd(w) )

    def AddConstraint( self,
                       np.ndarray[np.double_t, ndim=2, mode="c"] C,
                       np.ndarray[np.double_t, ndim=1, mode="c"] lbC,
                       np.ndarray[np.double_t, ndim=1, mode="c"] ubC,
                       hard_constraint = True):
        cdef bool flag = False
        if hard_constraint is True:
            self.thisptr.AddConstraint( NumpyToMatrixXd(C),
                                        NumpyToVectorXd(lbC),
                                        NumpyToVectorXd(ubC))
        else:
            self.thisptr.AddConstraint( NumpyToMatrixXd(C),
                                        NumpyToVectorXd(lbC),
                                        NumpyToVectorXd(ubC),
                                        flag)


    def Reset(self):
        self.thisptr.Reset()

    def SolveWBC(self):
        return self.thisptr.SolveWBC()

    def GetResult(self):
        cdef cWBCKits.VectorXd result = cWBCKits.VectorXd()
        self.thisptr.GetResult(result)
        return VectorXdToNumpy(result)

    def SetOptions(self, PyOptions options):
        return self.thisptr.SetOptions(deref(options.thisptr))
    
    def GetCostOfQP(self):
        return self.thisptr.GetCostOfQP()

cdef class RHP:
    cdef cWBCKits.RHP *thisptr
    # construct function in python
    def __cinit__(self):
        self.thisptr = new cWBCKits.RHP()

    # destruct function in python
    def __dealloc__(self):
        del self.thisptr

    # declare class methods    
    def AddTask(self,
                np.ndarray[np.double_t, ndim=2, mode="c"] A,
                np.ndarray[np.double_t, ndim=1, mode="c"] b,
                np.ndarray[np.double_t, ndim=1, mode="c"] w=None):
        if w is None:
            self.thisptr.AddTask( NumpyToMatrixXd(A),
                                  NumpyToVectorXd(b) )
        else:
            self.thisptr.AddTask( NumpyToMatrixXd(A),
                                  NumpyToVectorXd(b),
                                  NumpyToVectorXd(w) )

    def AddConstraint( self,
                       np.ndarray[np.double_t, ndim=2, mode="c"] C,
                       np.ndarray[np.double_t, ndim=1, mode="c"] lbC,
                       np.ndarray[np.double_t, ndim=1, mode="c"] ubC,
                       hard_constraint = True):
        cdef bool flag = False
        if hard_constraint is True:
            self.thisptr.AddConstraint( NumpyToMatrixXd(C),
                                        NumpyToVectorXd(lbC),
                                        NumpyToVectorXd(ubC))
        else:
            self.thisptr.AddConstraint( NumpyToMatrixXd(C),
                                        NumpyToVectorXd(lbC),
                                        NumpyToVectorXd(ubC),
                                        flag)


    def Reset(self):
        self.thisptr.Reset()

    def SetPriority(self,
                    np.ndarray[np.double_t, ndim=2, mode="c"] phi,
                    chi):
        cdef cWBCKits.VectorXi c_chi = cWBCKits.VectorXi(len(chi))
        for i in range(len(chi)):
            c_chi[i] = chi[i]

        self.thisptr.SetPriority( NumpyToMatrixXd(phi), c_chi)

    def SolveWBC(self):
        return self.thisptr.SolveWBC()

    def GetResult(self):
        cdef cWBCKits.VectorXd result = cWBCKits.VectorXd()
        self.thisptr.GetResult(result)
        return VectorXdToNumpy(result)

    def SetOptions(self, PyOptions options):
        return self.thisptr.SetOptions(deref(options.thisptr))
    
    def GetCostOfQP(self):
        return self.thisptr.GetCostOfQP()

cdef class HQP:
    cdef cWBCKits.HQP *thisptr
    # construct function in python
    def __cinit__(self):
        self.thisptr = new cWBCKits.HQP()

    # destruct function in python
    def __dealloc__(self):
        del self.thisptr

    # declare class methods    
    def AddTask(self,
                np.ndarray[np.double_t, ndim=2, mode="c"] A,
                np.ndarray[np.double_t, ndim=1, mode="c"] b,
                np.ndarray[np.double_t, ndim=1, mode="c"] w=None):
        if w is None:
            self.thisptr.AddTask( NumpyToMatrixXd(A),
                                  NumpyToVectorXd(b) )
        else:
            self.thisptr.AddTask( NumpyToMatrixXd(A),
                                  NumpyToVectorXd(b),
                                  NumpyToVectorXd(w) )

    def AddConstraint( self,
                       np.ndarray[np.double_t, ndim=2, mode="c"] C,
                       np.ndarray[np.double_t, ndim=1, mode="c"] lbC,
                       np.ndarray[np.double_t, ndim=1, mode="c"] ubC,
                       hard_constraint = True):
        cdef bool flag = False
        if hard_constraint is True:
            self.thisptr.AddConstraint( NumpyToMatrixXd(C),
                                        NumpyToVectorXd(lbC),
                                        NumpyToVectorXd(ubC))
        else:
            self.thisptr.AddConstraint( NumpyToMatrixXd(C),
                                        NumpyToVectorXd(lbC),
                                        NumpyToVectorXd(ubC),
                                        flag)


    def Reset(self):
        self.thisptr.Reset()

    def SetPriority(self,
                    phi,
                    chi):
        cdef cWBCKits.VectorXi c_phi = cWBCKits.VectorXi(len(phi))
        cdef cWBCKits.VectorXi c_chi = cWBCKits.VectorXi(len(chi))
        for i in range(len(phi)):
            c_phi[i] = phi[i]
        for i in range(len(chi)):
            c_chi[i] = chi[i]

        self.thisptr.SetPriority( c_phi, c_chi)

    def SolveWBC(self):
        return self.thisptr.SolveWBC()

    def GetResult(self):
        cdef cWBCKits.VectorXd result = cWBCKits.VectorXd()
        self.thisptr.GetResult(result)
        return VectorXdToNumpy(result)
    
    def SetOptions(self, PyOptions options):
        return self.thisptr.SetOptions(deref(options.thisptr))
    
    def GetCostOfQP(self):
        return self.thisptr.GetCostOfQP()