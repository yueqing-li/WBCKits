/**
 * @file wbc_base.h
 * @brief Abstract class of whole body control algorithm
 * @version 0.1
 * @date 2022-05 - 2022-09
 *  
 */
#ifndef WBCKITS_WBC_BASE_H_
#define WBCKITS_WBC_BASE_H_

#include <vector>
#include <iostream>
#include <qpOASES.hpp>
#include "task.h"
#include "constraint.h"
#include "options.h"


BEGIN_NAMESPACE_WBCKITS

class WBCBase
{
 public:
    WBCBase();

    virtual ~WBCBase();

    /**
     * @brief Solve the WBC problem using any kind of algorithm
     */
    virtual bool SolveWBC() = 0;

    /**
     * @brief add a subtask `min ||Ax-b||` to the task list
     */
    void AddTask(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);
    /**
     * @brief add a subtask `min ||w*(Ax-b)||` to the task list,
     *        here diagnore matrix w is stored as a vector
     */
    void AddTask(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
                 const Eigen::VectorXd &w);
    /**
     * @brief add a subconstrait `lbC <= C*x <= ubC` to the constraint list,
     * 
     * @param hard_constraint : default true, if set false, a slack variable will
     *                          be added to the constrait to form like 
     *                          `min ||w*v|| s.t. lbC <= C*x + v <= ubC`
     */
    void AddConstraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
                       const Eigen::VectorXd &ubC, bool hard_constraint = true);
   
    void AddConstraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
                       const Eigen::VectorXd &ubC, bool hard_constraint,
                       const Eigen::VectorXd &w);

    /**
     * @brief Reset WBC for solving a new WBC problem
     * 
     */
    void Reset();

    //set options with Options class
    virtual bool SetOptions( const Options& options);

    /**
     * @brief show tasks and constraints infomation
     *        TODO: combined with WBCOptions
     */
    virtual void ShowWBCInfo();

    inline int TasksSize() const { return task_list_.size(); }

    inline int ConstraintsSize() const { return constraint_list_.size(); }

    inline double GetCostOfQP() const {return qp_computation_cost_; }


 protected:
    std::vector<Task *> task_list_;
    std::vector<Constraint *> constraint_list_;
    Options options_;  // WBC options

    // ------------------ QP solver & settings ----------------------
    qpOASES::SQProblem* qp_ = NULL;
    qpOASES::Options qp_option_;
    qpOASES::returnValue status_code_solving_;
    bool init_done_;
    int nwsr_;
    double cpu_time_;
    int qp_var_;
    int qp_con_;

    bool new_qp_solver_;
    double qp_computation_cost_;  //> sum of time in solving qp
    Eigen::VectorXd opt_primal_last_;  //> optimal result from last solve    

    // ------------------------- Costs/Objects ----------------------
    /// Hessian matrix, H = A^T * A
    Eigen::MatrixXd hessian_mat_; 
    /// gradient vector, g = - A^T * b
    Eigen::VectorXd gradient_vec_; 
    /// weight
    Eigen::VectorXd omega_vec_;    
    /// A of final qp
    Eigen::MatrixXd task_mat_;     
    /// b of final qp
    Eigen::VectorXd target_vec_;
    //TODO config
    // epsilon regulation, if regulation used the final QP is 
    // `min ||Ax-b||+epsilon*||x||
    bool regulation_;
    double epsilon_;
    // Dimension
    int dim_var_{0};   ///< dimension of optimal variables = A.cols()
    int dim_obj_{0};   ///< dimension of all Tasks
    int dim_con_{0};   ///< dimension of all Constraints

    // ----------------- Bounds & Constraints -------------------------

    /// lb,ub TODO:not used yet
    double lbx_{-1.0e5};
    double ubx_{1.0e5};
    Eigen::VectorXd bound_vec_lbx_; 
    Eigen::VectorXd bound_vec_ubx_; 
    /// C
    Eigen::MatrixXd constraint_mat_; 
    Eigen::MatrixXd constraint_mat_t_;
    /// lbC,ubC
    Eigen::VectorXd constraint_vec_lbC_;
    Eigen::VectorXd constraint_vec_ubC_;

    void ClearTasks();
    void ClearConstraints();
    bool setQPOptions();
    bool createNewQPSolver();

};


END_NAMESPACE_WBCKITS




#endif  // WBC_BASE_H_