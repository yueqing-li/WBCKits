/**
 * @file task.h
 * @brief 
 * @version 0.1
 * @date 2022-05-11
 * 
 */

#ifndef WBCKITS_TASK_H_
#define WBCKITS_TASK_H_

#include "types.hpp"
#include <eigen3/Eigen/Dense>


BEGIN_NAMESPACE_WBCKITS
    
class Task
{
 public:
    Task();
    /**
     * @brief Construct a new Task object minimize ||Ax-b||
     */
    Task(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);
    /**
     * @brief Construct a new Task object minimize w||Ax-b||,
     *        here w as a diagnore matrix, is stored as a vector
     */
    Task(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
             const Eigen::VectorXd &w);

    ~Task() = default;
   

    inline Eigen::MatrixXd GetMat() const {return A_; }

    inline Eigen::VectorXd GetVec() const {return b_; }

    inline Eigen::VectorXd GetWet() const {return w_; }

    inline int TaskDim() const {return task_dim_; }

    inline int VarDim() const {return var_dim_; }


 protected:

    int task_dim_{0};         // the dimension of task, determining the row of A, b, w, ref
    int var_dim_{0};          // the DoF of variables in the WBC problem, determining the column of A
    // Task expression: min  W*( A*x - b ) //
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    Eigen::VectorXd w_;

};

END_NAMESPACE_WBCKITS


#endif  // TASK_BASE_H_