
#include "WBCKits/wbc_base.h"

BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/

WBCBase::WBCBase() : new_qp_solver_(true), init_done_(false), qp_var_(0), qp_con_(0)
{
    task_list_.clear();
    constraint_list_.clear();
}

WBCBase::~WBCBase()
{
    ClearTasks();
    ClearConstraints();
    delete qp_;
    qp_ = NULL;
}

void WBCBase::AddTask(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
{
    auto newTask = new Task(A, b);
    task_list_.push_back(newTask);

    dim_obj_ = dim_obj_ + newTask->TaskDim();
    if (dim_var_ == 0) {
        dim_var_ = newTask->VarDim();
    } else {
        assert(dim_var_ == newTask->VarDim());
    }
}

void WBCBase::AddTask(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
                      const Eigen::VectorXd &w)
{
    auto newTask = new Task(A, b, w);
    task_list_.push_back(newTask);

    dim_obj_ = dim_obj_ + newTask->TaskDim();
    if (dim_var_ == 0) {
        dim_var_ = newTask->VarDim();
    } else {
        assert(dim_var_ == newTask->VarDim());
    }
}

void WBCBase::AddConstraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
                            const Eigen::VectorXd &ubC, bool hard_constraint)
{
    constraint_list_.push_back(new Constraint(C, lbC, ubC, hard_constraint));
    dim_con_ = dim_con_ + C.rows();
}

void WBCBase::AddConstraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
                            const Eigen::VectorXd &ubC, bool hard_constraint,
                            const Eigen::VectorXd &w)
{
    constraint_list_.push_back(new Constraint(C, lbC, ubC, hard_constraint, w));
    dim_con_ = dim_con_ + C.rows();
}

bool WBCBase::SetOptions( const Options& options )
{
    options_ = options;

    return true;
}

void WBCBase::ShowWBCInfo()
{
    std::cout<<"Tasks' size is "<<task_list_.size()<<std::endl
             <<"Constraints' size is "<<constraint_list_.size()<<std::endl;
}

void WBCBase::Reset()
{
    ClearTasks();
    ClearConstraints();
}

/********************************
 *  P R O T E C T E D           *
 * ******************************/

void WBCBase::ClearTasks()
{
    for(auto taskPtr : task_list_) {
        delete taskPtr;
        taskPtr = NULL;
    }    
    task_list_.clear();
    dim_var_ = 0;
    dim_obj_ = 0;
}

void WBCBase::ClearConstraints()
{
    for(auto constrPtr : constraint_list_) {
        delete constrPtr;
        constrPtr = NULL;
    }
    constraint_list_.clear();
    dim_con_ = 0;
}

bool WBCBase::setQPOptions()
{
    if ( qp_ == NULL ) {
        return false;
    }

    qp_option_.printLevel = qpOASES::PL_LOW;
    // qp_option_.printLevel = qpOASES::PL_NONE;

    qp_->setOptions(qp_option_);

    return true;
}

bool WBCBase::createNewQPSolver()
{
    delete qp_;
    qp_ = new qpOASES::SQProblem(qp_var_, qp_con_, qpOASES::HST_SEMIDEF);

    if ( !setQPOptions() ) return false;

    init_done_ = false;
    new_qp_solver_ = false;

    return true;
}


END_NAMESPACE_WBCKITS

