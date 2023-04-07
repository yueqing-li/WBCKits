#include "WBCKits/wqp.h"

BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/

WQP::WQP()
{
    return;
}

WQP::~WQP()
{
    ClearTasks();
    ClearConstraints();
    delete qp_;
    qp_ = NULL;
}

bool WQP::SolveWBC()
{

    constructQP();
    if ( new_qp_solver_ )
        createNewQPSolver();
    
    return solveWQP();
}

void WQP::GetResult(Eigen::VectorXd &optimal_result)
{
    primal_opt_ = Eigen::VectorXd::Zero(qp_var_);

    qpOASES::real_t *primalPtr = new qpOASES::real_t[qp_var_];
    qp_->getPrimalSolution(primalPtr);

    for (int k=0; k<qp_var_; k++) {
        primal_opt_[k] = (double)primalPtr[k];
    }
    delete primalPtr;
    primalPtr = NULL;
    
    if (options_.get_slack_variable_) {
        optimal_result = primal_opt_;
    } else {
        optimal_result = primal_opt_.segment(0, dim_var_);
    }
}


/********************************
 *  P R I V A T E               *
 * ******************************/

bool WQP::constructQP()
{
    //select soft constraints and hard constraints
    Eigen::MatrixXd Cs, Ch;
    Eigen::VectorXd lds, uds, ws, ldh, udh;
    Cs.resize(dim_con_, dim_var_);
    lds.resize(dim_con_);
    uds.resize(dim_con_);
    ws.resize(dim_con_);
    Ch.resize(dim_con_, dim_var_);
    ldh.resize(dim_con_);
    udh.resize(dim_con_);
    int hard_c = 0;
    int soft_c = 0;
    for (auto constrPtr : constraint_list_) {
        int dim = constrPtr->ConstraintDim();
        if (constrPtr->IsHard()) {
            Ch.block(hard_c,0,dim,dim_var_) = constrPtr->GetMat();
            ldh.segment(hard_c,dim) = constrPtr->GetLower();
            udh.segment(hard_c,dim) = constrPtr->GetUpper();
            hard_c += dim;
        } else {
            Cs.block(soft_c,0,dim,dim_var_) = constrPtr->GetMat();
            lds.segment(soft_c,dim) = constrPtr->GetLower();
            uds.segment(soft_c,dim) = constrPtr->GetUpper();
            ws.segment(soft_c,dim) = constrPtr->GetWeight();
            soft_c += dim;
        }
    }

    Cs.conservativeResize(soft_c, dim_var_);
    lds.conservativeResize(soft_c);
    uds.conservativeResize(soft_c);
    ws.conservativeResize(soft_c);

    Ch.conservativeResize(hard_c, dim_var_);
    ldh.conservativeResize(hard_c);
    udh.conservativeResize(hard_c);

#ifdef __DEBUG__
    std::cout << "qp_var_ = " << qp_var_
              << ", (dim_var_ + soft_c) = " << (dim_var_ + soft_c) << std::endl
              << "qp_con_ = " << qp_con_
              << ", dim_con_ = " << dim_con_ << std::endl;
#endif 

    if ( qp_var_ != (dim_var_ + soft_c) ) {
        qp_var_ = dim_var_ + soft_c;
        new_qp_solver_ = true;
    }
    if ( qp_con_ != dim_con_ ) {
        qp_con_ = dim_con_;
        new_qp_solver_ = true;
    }
    // init A_all, b_all, w_all size, H,f
    Eigen::MatrixXd A_all;
    Eigen::VectorXd b_all, w_all;
    A_all   = Eigen::MatrixXd::Zero(dim_obj_, dim_var_);
    b_all = Eigen::VectorXd::Zero(dim_obj_);
    w_all  = Eigen::VectorXd::Zero(dim_obj_);
    
    // get A_all, b_all, w_all
    int start_row = 0;
    for (auto taskPtr : task_list_) {
        A_all.block(start_row,0, taskPtr->TaskDim(), dim_var_) = taskPtr->GetMat();
        b_all.segment(start_row, taskPtr->TaskDim()) = taskPtr->GetVec();
        w_all.segment(start_row, taskPtr->TaskDim()) = taskPtr->GetWet();

        start_row += taskPtr->TaskDim();
    }
    // times weight matrix
    A_all   = w_all.asDiagonal()*A_all;
    b_all = w_all.asDiagonal()*b_all;

    // calc H,f
    task_mat_.resize(dim_obj_+soft_c, qp_var_);
    task_mat_.setZero();
    task_mat_.block(0,0, dim_obj_, dim_var_) = A_all;
    task_mat_.block(dim_obj_, dim_var_, soft_c, soft_c) = ws.asDiagonal();
    target_vec_.resize(dim_obj_+soft_c);
    target_vec_.setZero();
    target_vec_.segment(0, dim_obj_) = b_all;

    hessian_mat_  = task_mat_.transpose()*task_mat_;
    gradient_vec_ = - task_mat_.transpose()*target_vec_;
    // regulation
    regulation_ = options_.enable_regulation_;
    epsilon_ = options_.weight_regulation_;
    if (regulation_) {
        hessian_mat_.block(0,0,dim_var_,dim_var_) =
        hessian_mat_.block(0,0,dim_var_,dim_var_) + epsilon_*Eigen::MatrixXd::Identity(dim_var_,dim_var_);
    }

    //init C_all, C_all^T, lbC, ubC
    constraint_mat_.resize(soft_c+hard_c, dim_var_+soft_c);
    constraint_mat_.setZero();
    constraint_mat_.block(0,0, soft_c, dim_var_) = Cs;
    constraint_mat_.block(0,dim_var_, soft_c, soft_c) = Eigen::MatrixXd::Identity(soft_c,soft_c);
    constraint_mat_.block(soft_c,0, hard_c, dim_var_) = Ch;

    constraint_vec_lbC_.resize(soft_c+hard_c);
    constraint_vec_lbC_.segment(0,soft_c)   = lds;
    constraint_vec_lbC_.segment(soft_c,hard_c) = ldh;

    constraint_vec_ubC_.resize(soft_c+hard_c);
    constraint_vec_ubC_.segment(0,soft_c)   = uds;
    constraint_vec_ubC_.segment(soft_c,hard_c) = udh;
/**
    constraint_mat_     = Eigen::MatrixXd::Zero(dim_con_, dim_var_);
    constraint_mat_t_   = Eigen::MatrixXd::Zero(dim_var_, dim_con_);
    constraint_vec_lbC_ = Eigen::VectorXd::Zero(dim_con_);
    constraint_vec_ubC_ = Eigen::VectorXd::Zero(dim_con_);
    //get C_all, lbC, ubC
    start_row = 0;
    for (auto constrPtr : constraint_list_) {
        constraint_mat_.block(start_row,0, constrPtr->ConstraintDim(), dim_var_)
            = constrPtr->GetMat();
        constraint_vec_lbC_.segment(start_row,constrPtr->ConstraintDim())
            = constrPtr->GetLower();
        constraint_vec_ubC_.segment(start_row,constrPtr->ConstraintDim())
            = constrPtr->GetUpper();

        start_row += constrPtr->ConstraintDim();
    }
**/
    constraint_mat_t_ = constraint_mat_.transpose();

    return true;
}

// bool WQP::createNewQPSolver()
// {
//     delete qp_;
//     qp_ = new qpOASES::SQProblem(qp_var_, dim_con_, qpOASES::HST_SEMIDEF);

//     if ( !setQPOptions() ) return false;

//     init_done_ = false;

//     return true;
// }

bool WQP::solveWQP()
{
    assert( qp_ != NULL );
    if (!init_done_) {
        nwsr_ = options_.nwsr_max_;
        cpu_time_ = options_.cpu_time_max_;
        status_code_solving_ = qp_->init(hessian_mat_.data(), gradient_vec_.data(),
                                         constraint_mat_t_.data(),
                                         bound_vec_lbx_.data(), bound_vec_ubx_.data(),
                                         constraint_vec_lbC_.data(), constraint_vec_ubC_.data(),
                                         nwsr_, &cpu_time_);

        if (status_code_solving_ > 0) {
            std::cout << "WQP solve failed!"<<".\n"
                  << "Error Message:"
                  << qpOASES::MessageHandling::getErrorCodeMessage(status_code_solving_)
                  << std::endl;
            qp_->reset();
            init_done_ = false;
            return false;
        } else {
            init_done_ = true;
        }
    } else {
#ifdef __DEBUG__
    printf("use hotstart in wqp.\n");
#endif
        nwsr_ = options_.nwsr_max_;
        cpu_time_ = options_.cpu_time_max_;
        status_code_solving_ = qp_->hotstart(hessian_mat_.data(), gradient_vec_.data(),
                                             constraint_mat_t_.data(),
                                             bound_vec_lbx_.data(), bound_vec_ubx_.data(),
                                             constraint_vec_lbC_.data(), constraint_vec_ubC_.data(),
                                             nwsr_, &cpu_time_);
        if (status_code_solving_ > 0) {
            std::cout << "WQP solve failed!"<<".\n"
                      << "Error Message:"
                      << qpOASES::MessageHandling::getErrorCodeMessage(status_code_solving_)
                      << std::endl;
            qp_->reset();
            init_done_ = false;
            return false;
        }
    }
    qp_computation_cost_ = cpu_time_;
    return true;
}


END_NAMESPACE_WBCKITS
