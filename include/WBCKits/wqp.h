/**
 * @file wqp.h
 * @brief implemention of weighted QP algorithm to solve WBC problem
 * @version 0.1
 * @date 2022-06 - 2022-09
 */

#ifndef WBCKITS_WQP_H_
#define WBCKITS_WQP_H_

#include "wbc_base.h"

BEGIN_NAMESPACE_WBCKITS

class WQP : public WBCBase
{
 public:
    WQP();

    virtual ~WQP();

    /**
     * @brief solve the WBC using WQP after added tasks and constraints
     * 
     * @return true : solve succeeded
     * @return false : solve failed
     */
    virtual bool SolveWBC() override;

    /**
     * @brief Get the Optimal result after solve the WQP
     * 
     * @param optimal_result : output, store the optimal result
     */
    void GetResult(Eigen::VectorXd & optimal_result);

 private:
    // construct the QP to solver from tasks list and constraints list
    bool constructQP();

   //  bool createNewQPSolver();

    bool solveWQP();

    Eigen::VectorXd primal_opt_;

};

END_NAMESPACE_WBCKITS


#endif  //WPQ_H_