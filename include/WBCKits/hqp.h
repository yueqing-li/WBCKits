/**
 * @file hqp.h
 * @brief implemention of HQP algorithm to solve WBC problem
 * @version 0.1
 * @date 2022-06-07
 * 
 */

#ifndef WBCKITS_HQP_H_
#define WBCKITS_HQP_H_

#include "wbc_base.h"

BEGIN_NAMESPACE_WBCKITS

class HQP : public WBCBase
{
 public:
    HQP();
    virtual ~HQP();

    /**
     * @brief Set the Priority order of tasks and constraints
     * 
     * @param phi priority vector of tasks, the order should be the same as that of AddTask()
     * @param chi priority vector of constraints, the order should be the same as that of AddConstraint()
     */
    void SetPriority(const Eigen::VectorXi &phi, const Eigen::VectorXi &chi);

    /**
     * @brief solve the WBC using HQP after adding tasks , constraints 
     *        and setting priority vector
     * @return true : solve succeeded
     * @return false : solve failed
     */
    virtual bool SolveWBC() override;

    /**
     * @brief Get the Optimal result after solve the WBC
     * 
     * @param optimal_result  output: store the optimal result
     */
    void GetResult(Eigen::VectorXd & optimal_result);

 protected:
// ---------------------  values -------------

    Eigen::VectorXi phi_;          // tasks priority matrix
    Eigen::VectorXi chi_;          // constraints priority vector
    Eigen::MatrixXd NullBase_;    // nullspace base matrix of higher level, shape: (dim_var_, remained DoF)

    Eigen::VectorXd xstar_;        // optimal result in level i
    // constraints from higher level
    Eigen::MatrixXd C_;
    Eigen::VectorXd C_ld_,C_ud_;
    // task at level i
    Eigen::MatrixXd Ai;
    Eigen::VectorXd bi, wi;
    // constraint at level i
    Eigen::MatrixXd Cis, Cih;
    Eigen::VectorXd lds, uds, wis, ldh, udh;

    int dim_u_; // variable dim without soft constraits
    int nv_i_;  // num of variable in level i, if constraint is hard, nv_i_ = dim_u_
    int nc_i_;  // num of constraint in level i

// ---------------------   functions  ------------------
    void initHQP();

    bool solveHQP();

    /**
     * @brief give sorted augmented task matrixs from Phi and tasks list in
     * priority i
     * @param i priority level, the highest level is starting from 1
     * @param Aia return value: augmented task matrixs
     * @param bia return value: augmented task target
     * @param wia return value: augmented weighted vector
     * @return false when no task select in this level
     */
    bool selectTask(int i, Eigen::MatrixXd &Aia, Eigen::VectorXd &bia,Eigen::VectorXd &wia);

    /*! 
     * @brief compute the orthonormal basis from augmented Jacobian
     * @param i  priority i
     * @param epsilon  threshold for judging rank
     * @return false when NullBasis is O
     */
    bool getNullBasis(int i, double eps = 1e-5);


    /**
     * @brief select constraints in priority level i
     * @param i priority level ,the highest level is starting from 1
     * @param Cis return : augmented soft constraint matrix in priority i
     * @param lds  return : lower bound of soft constraint
     * @param uds  return : upper bound of soft constraint
     * @param wis  return : weight of slack variable in soft constrains
     * @param Cih return : augmented hard constraint matrix in priority i
     * @param ldh  return : lower bound of hard constraint
     * @param udh  return : upper bound of hard constraint
     * @return false when no constraint in this level
     */
    bool selectConstraint(int i, Eigen::MatrixXd &Cis, Eigen::VectorXd &lds, Eigen::VectorXd &uds,
           Eigen::VectorXd &wis, Eigen::MatrixXd &Cih, Eigen::VectorXd &ldh, Eigen::VectorXd &udh);

    /**
     * @brief construct the QP in priority level i, 
     * @param i priority level ,the highest level is starting from 1
     * @return false when no task and constraint solve in this level
     */
    bool constructQP(int i);

    // solve the QP constructed from method ConstructQP()
    bool solveQP(int i);

    // bool createNewQPSolver();

    // debug use
    void showQPInfo(int i);
};

END_NAMESPACE_WBCKITS


#endif  // HQP_H_
