/**
 * @file constraint_base.h
 * @brief 
 * @version 0.1
 * @date 2022-05-11
 * 
 */

#ifndef WBCKITS_CONSTRAINT_H_
#define WBCKITS_CONSTRAINT_H_

#include "types.hpp"
#include <eigen3/Eigen/Dense>

BEGIN_NAMESPACE_WBCKITS
    
class Constraint
{
 public:
    Constraint();
    /**
     * @brief Construct a new Constraint object lbC <= C*x <= ubC
     *  
     * @param hard_constraint : default true, if set false, a slack variable will
     *                          be added to the constrait to form lbC <= C*x + v <= ubC
     */
    Constraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
               const Eigen::VectorXd &ubC, bool hard_constraint = true);
   
    Constraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
               const Eigen::VectorXd &ubC, bool hard_constraint,
               const Eigen::VectorXd &w);

    ~Constraint() = default;

   
    inline Eigen::MatrixXd GetMat() const {return C_; }

    inline Eigen::VectorXd GetLower() const {return lbC_; }

    inline Eigen::VectorXd GetUpper() const {return ubC_; }

    inline Eigen::VectorXd GetWeight() const {return w_; }

    inline int ConstraintDim() const {return constr_dim_; }

    inline int VarDim() const {return var_dim_; }

    inline bool IsHard() const {return hard_constraint_; }


 protected:

    int constr_dim_{0};             // the dimension of constraint, determining the row of C, lbC, ubC
    int var_dim_{0};         // the DoF of variables in the optimization, determining the column of C
    // Constraint expression: lbC <= C*x <= ubC //
    Eigen::MatrixXd C_;
    Eigen::VectorXd lbC_;
    Eigen::VectorXd ubC_;

    Eigen::VectorXd w_;  // weight vector for solf constrait lbC <= C*x + v <= ubC 

    bool hard_constraint_ = true;

};

END_NAMESPACE_WBCKITS



#endif // CONSTRAINT_BASE_H_
