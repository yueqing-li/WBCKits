/**
 * @file options.h
 * @brief settings for all wbc solver class
 * @version 0.1
 * @date 2022-08-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef WBCKITS_OPTIONS_H_
#define WBCKITS_OPTIONS_H_

#include "types.hpp"

BEGIN_NAMESPACE_WBCKITS

class Options
{
// public functions
 public:
    
    Options();
    // deep copy
    Options( const Options & rhs );

    ~Options();

    // deep copy
    Options& operator= ( const Options& rhs);

    bool SetDefault();

    void Print() const;

 protected:
    bool copy( const Options& rhs);

// public values
 public:
    
    // ------------------ QP solver settings ----------------------
    int nwsr_max_;  // Maximum number of working set recalculations
    double cpu_time_max_;  // max time used for solving a single qp

    // -------------------- QP Construction ----------------------
    bool enable_regulation_;
    double weight_regulation_;

   // -------------------- HQP  ----------------------
    AlgorithmOfHQP hqp_algorithm_;
    DecomposeMethod decompose_method_;
    double decompose_threshold_;

    // -------------------- Optimal Results -----------------------
    bool get_slack_variable_;


};


END_NAMESPACE_WBCKITS


#endif  // WBCKITS_OPTIONS_H_