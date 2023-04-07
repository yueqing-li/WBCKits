/**
 * @file rhp.h
 * @brief implemention of RHP-HQP algorithm to solve WBC problem
 * @version 0.1
 * @date 2022-06-07
 * 
 */

#ifndef WBCKITS_RHP_H_
#define WBCKITS_RHP_H_

#include "wbc_base.h"

BEGIN_NAMESPACE_WBCKITS

class RHP : public WBCBase
{
 public:
    RHP();

    virtual ~RHP();

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

    Eigen::MatrixXd phi_;          // tasks priority matrix
    Eigen::VectorXi chi_;          // constraints priority vector
    Eigen::MatrixXd Projector_;    // projector to higher level

    Eigen::VectorXd xstar_;        // optimal result in level i
    // constraints from higher level
    Eigen::MatrixXd C_;
    Eigen::VectorXd C_ld_,C_ud_;

    int nv_i_;  // num of variable in level i, if constraint is hard, nv_i_ = dim_var_
    int nc_i_;  // num of constraint in level i

    
 // ---------------------   functions  ------------------
    void initRHP();

    bool solveRHP();

    /**
     * @brief give sorted augmented task matrixs from Phi and tasks list in
     * priority i
     * @param i priority level, the highest level is starting from 1
     * @param Aisa return value: sorted augmented task matrixs
     * @param bisa return value: sorted augmented task target
     * @param wisa return value: sorted augmented weighted vector
     * @param Lambdai return value: phi(i,s)
     * @return false when no task select in this level
     */
    bool selectAndSortTask(int i, Eigen::MatrixXd &Aisa, Eigen::VectorXd &bisa,Eigen::VectorXd &wisa, Eigen::VectorXd &Lambdai, double eps = 1e-5);

    /*! 
     * @brief compute the orthonormal basis from augmented Jacobian
     * @param Aisa  augmented Jacobian sorted by matrix Ai
     * @param B  output : orthornormal basis
     * @param origin  output : basis index from J
     * @param r  output : the raw number of B
     * @param epsilon  threshold for avoiding norm of basis diveding zero
     */
    void getOrthBasis(const Eigen::MatrixXd &Aisa, Eigen::MatrixXd &B, Eigen::VectorXi &origin, int &r, double eps = 1e-5);

    /*! compute the generalized projector of priority level i
     * @param i  priority level, the highest level is starting from 0
     */
    void getGeneralizedProjector(int i);

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

   //  bool createNewQPSolver();

    // debug use
    void showQPInfo(int i);


};

END_NAMESPACE_WBCKITS



#endif  //RHP_H_