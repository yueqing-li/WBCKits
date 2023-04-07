

#ifndef WBCKITS_TYPES_HPP_
#define WBCKITS_TYPES_HPP_

/* Uncomment the following line to enable debug information. */
// #define __DEBUG__

#define BEGIN_NAMESPACE_WBCKITS namespace WBCKits {

#define END_NAMESPACE_WBCKITS   }

#define WBCKITS_ERROR false


BEGIN_NAMESPACE_WBCKITS

//> HQP construct form
enum AlgorithmOfHQP
{
    HQP_Origin,                 // using original form to solve HQP problem
    HQP_Nullspace               // using parameterized form with nullspace base to solve HQP
};

//> different decompose methods to calculate nullspace base
enum DecomposeMethod
{
    Decompose_QR = 0,                         // TODO
    Decompose_SVD,                        // using SVD
    Decompose_COD                         // TODO
};
    
END_NAMESPACE_WBCKITS


#endif // WBCKITS_TYPES_HPP_