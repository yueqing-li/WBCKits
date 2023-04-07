
#include "WBCKits/task.h"
#include <cassert>

BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/
Task::Task()
{
    //TODO: if add reset method, implement default function
    assert(WBCKITS_ERROR && "no task define here ");
}

Task::Task(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
{
    task_dim_ = A.rows();
    var_dim_ = A.cols();
    assert(task_dim_==b.size());

    A_ = A;
    b_ = b;

    w_ = Eigen::VectorXd::Ones(task_dim_);
}

Task::Task(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
           const Eigen::VectorXd &w)
{
    task_dim_ = A.rows();
    var_dim_ = A.cols();
    assert(task_dim_==b.size());
    assert(task_dim_==w.size());

    A_ = A;
    b_ = b;
    w_ = w;
}


END_NAMESPACE_WBCKITS



