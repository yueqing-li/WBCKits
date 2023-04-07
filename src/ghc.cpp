#include "WBCKits/ghc.h"

BEGIN_NAMESPACE_WBCKITS
/********************************
 *  P U B L I C                 *
 * ******************************/
GHC::GHC()
{
    return;
}

GHC::~GHC()
{
    ClearTasks();
    ClearConstraints();
    delete qp_;
    qp_ = NULL;
}

void GHC::SetPriority(const Eigen::MatrixXd &phi, const Eigen::VectorXi &chi)
{
    if (task_list_.size() != phi.cols()) {
        assert(WBCKITS_ERROR && "tasks num doesn't match list's size.");
    }
    if (constraint_list_.size() != chi.size()) {
        assert(WBCKITS_ERROR &&  "constraints num doesn't match list's size.");
    }

    phi_ = phi;
    chi_ = chi;
}

void GHC::SetPriority(const Eigen::VectorXi &phi, const Eigen::VectorXi &chi)
{
    if (task_list_.size() != phi.size()) {
        assert(WBCKITS_ERROR && "tasks num doesn't match list's size.");
    }
    if (constraint_list_.size() != chi.size()) {
        assert(WBCKITS_ERROR &&  "constraints num doesn't match list's size.");
    }

    Eigen::MatrixXd phi_mat;
    int nt = phi.size();
    phi_mat.resize(nt,nt);
    phi_mat.setZero();

    for ( int i=0; i<nt; i++ ) {
        for ( int j=0; j<nt; j++ ) {
            if ( (phi(i) > phi(j) && phi(j) != 0) || phi(i) == 0 )
            {
                phi_mat(i, j) = 1.0;  // task i priority lower than task j or task i has been removed
            } else {
                phi_mat(i, j) = 0.0;
            }
        }
    }

    phi_ = phi_mat;
    chi_ = chi;
}

bool GHC::SolveWBC()
{
    constructQP();
    createNewQPSolver();
    return solveGHC();
}

void GHC::GetResult(Eigen::VectorXd & optimal_result)
{
    primal_opt_ = Eigen::VectorXd::Zero(qp_var_);

    qpOASES::real_t *primalPtr = new qpOASES::real_t[qp_var_];
    qp_->getPrimalSolution(primalPtr);

    for (int k=0; k<qp_var_; k++) {
        primal_opt_[k] = (double)primalPtr[k];
    }
    delete primalPtr;
    primalPtr = NULL;

    augment_x_ = primal_opt_.segment(0, dim_var_*task_list_.size());
    optimal_result = GeneralP_*augment_x_;
}

/********************************
 *  P R I V A T E               *
 * ******************************/

bool GHC::constructQP()
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

    // A_all is diagonal block matrix diag(A1,A2,...,An)
    Eigen::MatrixXd A_all;
    Eigen::VectorXd b_all, w_all;
    int num_task = task_list_.size();
    A_all   = Eigen::MatrixXd::Zero(dim_obj_, dim_var_*num_task);
    b_all = Eigen::VectorXd::Zero(dim_obj_);
    // w_all  = Eigen::VectorXd::Zero(dim_obj_);  // task weight is not used in GHC
    augment_J_ = Eigen::MatrixXd::Zero(dim_obj_, dim_var_);
 
    // get A_all, b_all, w_all
    int start_row = 0;
    int start_col = 0;
    for (auto taskPtr : task_list_) {
        A_all.block(start_row, start_col, taskPtr->TaskDim(), dim_var_) = taskPtr->GetMat();
        b_all.segment(start_row, taskPtr->TaskDim()) = taskPtr->GetVec();
        // w_all.segment(start_row, taskPtr->TaskDim()) = taskPtr->GetWet();
        augment_J_.block(start_row, 0, taskPtr->TaskDim(), dim_var_) = taskPtr->GetMat();

        start_row += taskPtr->TaskDim();
        start_col += dim_var_;
    }

    // get H, f
    qp_var_ = dim_var_*num_task + soft_c;
    task_mat_.resize(dim_obj_+soft_c, qp_var_);
    task_mat_.setZero();
    task_mat_.block(0,0, dim_obj_, dim_var_*num_task) = A_all;
    task_mat_.block(dim_obj_, dim_var_*num_task, soft_c, soft_c) = ws.asDiagonal();
    target_vec_.resize(dim_obj_+soft_c);
    target_vec_.setZero();
    target_vec_.segment(0, dim_obj_) = b_all;

    hessian_mat_  = task_mat_.transpose()*task_mat_;
    gradient_vec_ = - task_mat_.transpose()*target_vec_;
    // regulation
    regulation_ = options_.enable_regulation_;
    epsilon_ = options_.weight_regulation_;
    if (regulation_) {
        hessian_mat_.block(0,0,dim_var_*num_task,dim_var_*num_task) 
            = hessian_mat_.block(0,0,dim_var_*num_task,dim_var_*num_task)
            + epsilon_*Eigen::MatrixXd::Identity(dim_var_*num_task,dim_var_*num_task);
    }

    // calculate Generalized Projector
    GeneralP_ = Eigen::MatrixXd::Zero(dim_var_, dim_var_*num_task);
    start_col = 0;
    Eigen::MatrixXd tempPi;
    for ( int i=0; i<num_task; i++) {
        getGeneralizedProjector(i, tempPi);
        GeneralP_.block(0, start_col, dim_var_, dim_var_) = tempPi;
        start_col += dim_var_;
    }

    //init constraint_mat_, constraint_mat_^T, lbC, ubC
    constraint_mat_.resize(soft_c+hard_c, qp_var_);
    constraint_mat_.setZero();
    constraint_mat_.block(0,0, soft_c, dim_var_*num_task) = Cs*GeneralP_;
    constraint_mat_.block(0,dim_var_*num_task, soft_c, soft_c) = Eigen::MatrixXd::Identity(soft_c,soft_c);
    constraint_mat_.block(soft_c,0, hard_c, dim_var_*num_task) = Ch*GeneralP_;

    constraint_vec_lbC_.resize(soft_c+hard_c);
    constraint_vec_lbC_.segment(0,soft_c)   = lds;
    constraint_vec_lbC_.segment(soft_c,hard_c) = ldh;

    constraint_vec_ubC_.resize(soft_c+hard_c);
    constraint_vec_ubC_.segment(0,soft_c)   = uds;
    constraint_vec_ubC_.segment(soft_c,hard_c) = udh;

    constraint_mat_t_ = constraint_mat_.transpose();

    return true;

}

bool GHC::getGeneralizedProjector(int i, Eigen::MatrixXd & Pi, double eps)
{
    int num_task = task_list_.size();
    Eigen::VectorXd alpha_i, alpha_i_sort;
    alpha_i.resize(dim_obj_);
    alpha_i_sort.resize(dim_obj_);
    int start_j = 0;
    for ( int kk=0; kk<num_task; kk++ ) {
        alpha_i.segment(start_j, task_list_[kk]->TaskDim()) = phi_(i, kk)*Eigen::VectorXd::Ones(task_list_[kk]->TaskDim());
        start_j += task_list_[kk]->TaskDim();
    }

    Eigen::VectorXi index;
    index.resize(dim_obj_);

    // init index
    for ( int kk=0; kk<dim_obj_; kk++ ) index(kk) = kk;
    // descending order of phi_.row(i)
    for ( int kk=0; kk<dim_obj_; kk++ ) {
        for ( int j=0; j<dim_obj_-kk-1; j++ ) {
            if ( alpha_i(j) < alpha_i(j+1) ) {
                int tmp = index(j);
                index(j) = index(j+1);
                index(j+1) = tmp;
            }
        }
    }
    for ( int kk=0; kk<dim_obj_; kk++ ) {
        alpha_i_sort(kk) = alpha_i(index(kk));
    }

    // get sorted augmented task matrix from alpha_i descending order
    Eigen::MatrixXd A_sort;
    A_sort.resize(dim_obj_, dim_var_);
    A_sort.setZero();
    for ( int kk=0; kk<dim_obj_; kk++ ) {
        A_sort.row(kk) = augment_J_.row(index(kk));
    }

    Eigen::MatrixXd OrthBasis;
    Eigen::VectorXd alpha_r;
    OrthBasis.resize(dim_obj_,dim_var_);
    alpha_r.resize(dim_obj_);
    int rr = 0;
    for ( int kk=0; kk<dim_obj_; kk++ ) {
        if ( rr >= dim_var_ ) break; // full colum rank, no more calc need

        OrthBasis.row(rr) = A_sort.row(kk);
        for ( int ss=0; ss<rr; ss++ ) {
            // orthogonalization
            OrthBasis.row(rr) = OrthBasis.row(rr) 
                        - (OrthBasis.row(rr)*(OrthBasis.row(ss).transpose()))*OrthBasis.row(ss);
        }
        double scalar = std::sqrt(OrthBasis.row(rr)*(OrthBasis.row(rr).transpose()));
        if ( scalar > eps ) {
            OrthBasis.row(rr) = OrthBasis.row(rr)/scalar;
            alpha_r(rr) = alpha_i_sort(kk);
            rr++;
        }
    }

    OrthBasis.conservativeResize(rr, dim_var_);
    alpha_r.conservativeResize(rr);

    // get Pi
    Pi.resize(dim_var_, dim_var_);
    Pi = Eigen::MatrixXd::Identity(dim_var_, dim_var_) 
                        - OrthBasis.transpose()*alpha_r.asDiagonal()*OrthBasis;

    return true;
}

bool GHC::createNewQPSolver()
{
    delete qp_;
    qp_ = new qpOASES::SQProblem(qp_var_, dim_con_, qpOASES::HST_SEMIDEF);

    if ( !setQPOptions() ) return false;

    init_done_ = false;

    return true;
}

bool GHC::solveGHC()
{
    if (!init_done_) {
        nwsr_ = options_.nwsr_max_;
        cpu_time_ = options_.cpu_time_max_;
        status_code_solving_ = qp_->init(hessian_mat_.data(), gradient_vec_.data(),
                                         constraint_mat_t_.data(),
                                         bound_vec_lbx_.data(), bound_vec_ubx_.data(),
                                         constraint_vec_lbC_.data(), constraint_vec_ubC_.data(),
                                         nwsr_, &cpu_time_);

        if (status_code_solving_ > 0) {
            std::cout << "GHC solve failed!"<<".\n"
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
        nwsr_ = options_.nwsr_max_;
        cpu_time_ = options_.cpu_time_max_;
        status_code_solving_ = qp_->hotstart(hessian_mat_.data(), gradient_vec_.data(),
                                             constraint_mat_t_.data(),
                                             bound_vec_lbx_.data(), bound_vec_ubx_.data(),
                                             constraint_vec_lbC_.data(), constraint_vec_ubC_.data(),
                                             nwsr_, &cpu_time_);
        if (status_code_solving_ > 0) {
            std::cout << "GHC solve failed!"<<".\n"
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

void GHC::showQPInfo()
{
    std::cout << "hessian:\n" << hessian_mat_ << std::endl
              << "gradient_vec :" << gradient_vec_.transpose() << std::endl
              << "target mat :\n" << task_mat_ << std::endl
              << "target vec : " << target_vec_.transpose() << std::endl;

    std::cout << "------------------------------------------------------" << std::endl;

    std::cout << "constraint_mat_:\n" << constraint_mat_t_.transpose() << std::endl
              << "constraint_vec_lbC_: " << constraint_vec_lbC_.transpose() << std::endl
              << "constraint_vec_ubC_: " << constraint_vec_ubC_.transpose() << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "generalization project matrix:" << std::endl
              << GeneralP_ << std::endl;
    std::cout << "--------------------- END ----------------------------" << std::endl;
    
}

END_NAMESPACE_WBCKITS