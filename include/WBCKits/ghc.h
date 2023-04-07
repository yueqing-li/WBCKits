/**
 * @file ghc.h
 * @brief implemention of GHC algorithm to solve WBC problem
 * @version 0.1
 * @date 2022-08 - 2022-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef WBCKITS_GHC_H_
#define WBCKITS_GHC_H_

#include "wbc_base.h"

BEGIN_NAMESPACE_WBCKITS

class GHC : public WBCBase
{
 public:
    GHC();
    virtual ~GHC();

    /**
     * @brief Set the Priority of tasks and constraints
     * 
     * @param phi priority matrix of tasks, the column order should be the same as that of AddTask()
     * @param chi priority vector of constraints, the order should be the same as that of AddConstraint()
     */
    void SetPriority(const Eigen::MatrixXd &phi, const Eigen::VectorXi &chi);

    /**
     * @brief Set the Priority of tasks and constraints similar to HQP
     * 
     * @param phi priority vector of tasks, the order should be the same as that of AddTask()
     * @param chi priority vector of constraints, the order should be the same as that of AddConstraint()
     */
    void SetPriority(const Eigen::VectorXi &phi, const Eigen::VectorXi &chi);

    /**
     * @brief solve the WBC using RHP after adding tasks , constraints 
     *        and setting priority matrix
     * @return true : solve succeeded
     * @return false : solve failed
     */
    virtual bool SolveWBC() override;

    /**
     * @brief Get the Optimal result after solve the WBC
     * 
     * @param optimal_result  output, store the optimal result
     */
    void GetResult(Eigen::VectorXd & optimal_result);

 protected:
    // ---------------------  values -------------

    Eigen::MatrixXd phi_;          //< tasks priority matrix
    Eigen::VectorXi chi_;          //< constraints priority vector
    Eigen::MatrixXd augment_J_;    //< augmented task matrix
    Eigen::MatrixXd GeneralP_;     //< augmented general projectors defined by phi_ 
    Eigen::VectorXd augment_x_;    //< augmented Optimization variables

    Eigen::VectorXd primal_opt_;

    // ---------------------  functions -------------
    bool constructQP();

    /**
     * @brief Get the Generalized Projector object after setting 
     * 
     * @param i  the i-th task
     * @param Pi output: the generalized projector of task i
     * @return false 
     */
    bool getGeneralizedProjector(int i, Eigen::MatrixXd & Pi, double eps = 1e-5);

    bool createNewQPSolver();

    bool solveGHC();
   
    void showQPInfo();

};


END_NAMESPACE_WBCKITS


#endif  // WBCKITS_GHC_H_