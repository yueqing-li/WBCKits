#include "WBCKits/hqp.h"


BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/

HQP::HQP()
{
    options_.enable_regulation_ = true;
}

HQP::~HQP()
{
    ClearTasks();
    ClearConstraints();
    delete qp_;
    qp_ = NULL;
}

void HQP::SetPriority(const Eigen::VectorXi &phi, const Eigen::VectorXi &chi)
{

    if (task_list_.size() != phi.size()) {
        assert(WBCKITS_ERROR && "tasks num doesn't match list's size.");
    }
    if (constraint_list_.size() != chi.size()) {
        assert(WBCKITS_ERROR && "constraints num doesn't match list's size." );
    }

    phi_ = phi;
    chi_ = chi;

}

bool HQP::SolveWBC()
{
    return solveHQP();
}

void HQP::GetResult(Eigen::VectorXd & optimal_result)
{
    optimal_result = xstar_;
}


/********************************
 *  P R I V A T E               *
 * ******************************/
void HQP::initHQP()
{
    xstar_ = Eigen::VectorXd::Zero(dim_var_);
    C_.resize(0, dim_var_);
    C_ld_.resize(0);
    C_ud_.resize(0);
    qp_computation_cost_ = 0.;
}

bool HQP::solveHQP()
{
    bool solve_flag = true;
    initHQP();
    int priority_level = phi_.maxCoeff();

    if ( options_.hqp_algorithm_ == HQP_Origin ) {  // solve HQP by Original form
        for ( int i=1; i<=priority_level; i++ ) {
            if ( constructQP(i) ) {
                if ( solveQP(i) ) {
                    qpOASES::real_t *primal_ptr = new qpOASES::real_t[nv_i_];
                    Eigen::VectorXd vstar, primal_opt;
                    primal_opt.resize(nv_i_);
                    qp_->getPrimalSolution(primal_ptr);
                    for (int k=0; k<nv_i_; k++) {
                        primal_opt(k) = (double)(primal_ptr[k]);
                    }
                    delete primal_ptr;

                    int nt_i = Ai.rows();
                    int ncs_i_new = nv_i_ - dim_u_;
                    xstar_ = primal_opt.segment(0, dim_u_);
                    vstar = primal_opt.segment(dim_u_, ncs_i_new);
                    
                    int nC = C_.rows();
                    if (nt_i != 0) {
                        C_ld_.segment(nC-ncs_i_new-nt_i, nt_i) = Ai*xstar_;
                        C_ud_.segment(nC-ncs_i_new-nt_i, nt_i) = Ai*xstar_;
                    }
                    if (ncs_i_new != 0) {
                        C_ld_.segment(nC-ncs_i_new, ncs_i_new) 
                                = C_ld_.segment(nC-ncs_i_new, ncs_i_new) - vstar;
                        C_ud_.segment(nC-ncs_i_new, ncs_i_new)
                                = C_ud_.segment(nC-ncs_i_new, ncs_i_new) - vstar;
                    }
                    
                } else {
                    solve_flag = false;
                }
            }
        }
    } else if ( options_.hqp_algorithm_ == HQP_Nullspace ) {  // solve HQP by parameterize form
        for ( int i=1; i<= priority_level; i++ ) {
            bool haveDoF =  getNullBasis(i-1, options_.decompose_threshold_);
            if (!haveDoF) { // no DoF left, algorithm finish
                return solve_flag;
            }

            if ( constructQP(i) ) {
                if ( solveQP(i) ) {
                    qpOASES::real_t *primal_ptr = new qpOASES::real_t[nv_i_];
                    Eigen::VectorXd ustar, vstar, primal_opt;
                    primal_opt.resize(nv_i_);
                    qp_->getPrimalSolution(primal_ptr);
                    for (int k=0; k<nv_i_; k++) {
                        primal_opt(k) = (double)(primal_ptr[k]);
                    }
                    delete primal_ptr;

                    int ncs_i_new = nv_i_ - dim_u_;
                    ustar = primal_opt.segment(0, dim_u_);
                    vstar = primal_opt.segment(dim_u_, ncs_i_new);

                    xstar_ = NullBase_*ustar + xstar_;
                    
                    if (ncs_i_new != 0) {
                        int nC = C_.rows();
                        C_ld_.segment(nC-ncs_i_new, ncs_i_new) 
                                = C_ld_.segment(nC-ncs_i_new, ncs_i_new) - vstar;
                        C_ud_.segment(nC-ncs_i_new, ncs_i_new)
                                = C_ud_.segment(nC-ncs_i_new, ncs_i_new) - vstar;
                    }
                    
                } else {
                    solve_flag = false;
                }
            }
        }
    }

    return solve_flag;

}

bool HQP::selectTask(int i, Eigen::MatrixXd &Aia, Eigen::VectorXd &bia,Eigen::VectorXd &wia)
{
    Aia.resize(dim_obj_, dim_var_);
    bia.resize(dim_obj_);
    wia.resize(dim_obj_);

    int nt = phi_.size();  // number of tasks
    int start_rows = 0;
    for ( int k=0; k < nt; k++ ) {
        if ( i == phi_[k] ) {
            int dim = task_list_[k]->TaskDim();
            Aia.block(start_rows, 0, dim, dim_var_) = task_list_[k]->GetMat();
            bia.segment(start_rows, dim) = task_list_[k]->GetVec();
            wia.segment(start_rows, dim) = task_list_[k]->GetWet();
            start_rows += dim;
        }
    }
    Aia.conservativeResize(start_rows, dim_var_);
    bia.conservativeResize(start_rows,1);
    wia.conservativeResize(start_rows,1);

    if ( start_rows == 0 ) return false;

    return true;
}

bool HQP::getNullBasis(int i, double eps)
{
    if ( i == 0 ) {
        NullBase_ = Eigen::MatrixXd::Identity(dim_var_, dim_var_);
        return true;
    }
    /**
     * TODO: decrease use of selectTask to get Aia
     */
    Eigen::MatrixXd Ai_hat, Aia, null_Ai_h;
    Eigen::VectorXd bi, wi;
    bool have_task;
    have_task = selectTask(i, Aia, bi, wi);
    if ( !have_task ) {
        return true; // no change to nullbase in this level
    }
    Ai_hat = Aia*NullBase_;
    if ( options_.decompose_method_ == Decompose_SVD ) {
        Eigen::BDCSVD<Eigen::MatrixXd> svd;
        svd.setThreshold(eps);  // singular < eps*max_singular -> 0
        svd.compute(Ai_hat, Eigen::ComputeThinU|Eigen::ComputeFullV);
        int rank_Ai_h = svd.rank();
        if ( Ai_hat.cols() == rank_Ai_h ) {  // Ai_hat full rank, nullBase is O
            return false;
        }
        null_Ai_h = svd.matrixV().rightCols( Ai_hat.cols() - svd.rank() );
        // std::cout << "SVD null A "<<i<<"_h\n" << null_Ai_h << std::endl;
    } else if ( options_.decompose_method_ == Decompose_QR ) {
        //TODO
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
        qr.setThreshold(eps);
        qr.compute(Ai_hat.transpose()*Ai_hat);
        int rank_Ai_h = qr.rank();
        if ( Ai_hat.cols() == rank_Ai_h ) {  // Ai_hat full rank, nullBase is O
            return false;
        }
        null_Ai_h = qr.householderQ();
        null_Ai_h = null_Ai_h.rightCols(Ai_hat.cols() - qr.rank() );
        // std::cout << "QR null Ai h\n" << null_Ai_h << std::endl;
        // Eigen::MatrixXd per;
        // per = qr.colsPermutation();
        // std::cout << "permutation matrix:\n" << per << std::endl;
    } else if ( options_.decompose_method_ == Decompose_COD ) {
        //TODO
        return false;
    }
    NullBase_ = NullBase_*null_Ai_h;

    return true;
}

bool HQP::selectConstraint(int i, Eigen::MatrixXd &Cis, Eigen::VectorXd &lds, Eigen::VectorXd &uds,
                           Eigen::VectorXd &wis, Eigen::MatrixXd &Cih, Eigen::VectorXd &ldh, Eigen::VectorXd &udh)
{
    Cis.resize(dim_con_, dim_var_);
    lds.resize(dim_con_);
    uds.resize(dim_con_);
    wis.resize(dim_con_);
    Cih.resize(dim_con_, dim_var_);
    ldh.resize(dim_con_);
    udh.resize(dim_con_);

    int nc = chi_.size();  // num of constraits
    int hard_c = 0;
    int soft_c = 0;
    for (int k=0; k<nc; k++) {
        if (i == chi_[k]) {
            int dim = constraint_list_[k]->ConstraintDim();
            if (constraint_list_[k]->IsHard()) {
                Cih.block(hard_c,0,dim,dim_var_) = constraint_list_[k]->GetMat();
                ldh.segment(hard_c,dim) = constraint_list_[k]->GetLower();
                udh.segment(hard_c,dim) = constraint_list_[k]->GetUpper();
                hard_c += dim;
            } else {
                Cis.block(soft_c,0,dim,dim_var_) = constraint_list_[k]->GetMat();
                lds.segment(soft_c,dim) = constraint_list_[k]->GetLower();
                uds.segment(soft_c,dim) = constraint_list_[k]->GetUpper();
                wis.segment(soft_c,dim) = constraint_list_[k]->GetWeight();
                soft_c += dim;
            }
        }
    }

    Cis.conservativeResize(soft_c, dim_var_);
    lds.conservativeResize(soft_c);
    uds.conservativeResize(soft_c);
    wis.conservativeResize(soft_c);

    Cih.conservativeResize(hard_c, dim_var_);
    ldh.conservativeResize(hard_c);
    udh.conservativeResize(hard_c);

    if ((soft_c+hard_c) == 0) return false;

    return true;
}

bool HQP::constructQP(int i)
{
    bool newQP = false;

    bool haveTask;
    haveTask = selectTask(i, Ai, bi, wi);

    bool haveConstraint;  
    haveConstraint  = selectConstraint(i, Cis, lds, uds, wis, Cih, ldh, udh);

    int ncs = 0;  // num of soft constraint
    int nch = 0;  // num of hard constraint
    int nt = 0;
    dim_u_ = 0;
    nv_i_ = 0;
    nc_i_ = 0;
    if (haveTask && haveConstraint) {
        ncs = Cis.rows();
        nch = Cih.rows();
        nt = Ai.rows();
    } else if (haveTask && !haveConstraint) {
        nt = Ai.rows();
    } else if (!haveTask && haveConstraint) {
        ncs = Cis.rows();
        nch = Cih.rows();
    } else if (!haveTask && !haveConstraint) {
        return false;
    }


    if ( options_.hqp_algorithm_ == HQP_Origin ) {  // use original form construct QP
        dim_u_ = dim_var_;

        if (nv_i_ != (dim_var_ + ncs)){
            nv_i_ = dim_var_ + ncs;
            newQP = true;
        }
        // task mat
        task_mat_.resize(nt+ncs, nv_i_);
        task_mat_.setZero();
        task_mat_.block(0,0, nt, dim_var_) = wi.asDiagonal()*Ai;
        // task_mat_.block(nt, dim_var_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
        task_mat_.block(nt, dim_var_, ncs, ncs) = wis.asDiagonal();

        target_vec_.resize(nt+ncs);
        target_vec_.setZero();
        target_vec_.segment(0, nt) =  wi.asDiagonal()*bi;

        int nc_pre = C_.rows();
        if (nc_i_ != (ncs + nch + nc_pre)) {
            nc_i_ = ncs + nch + nc_pre;
            newQP = true;
        }
        // constraint
        if (nc_pre == 0) {
            constraint_mat_.resize(ncs+nch, nv_i_);
            constraint_mat_.setZero();
            constraint_mat_.block(0,0, ncs, dim_var_) = Cis;
            constraint_mat_.block(0,dim_var_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
            constraint_mat_.block(ncs,0, nch, dim_var_) = Cih;

            constraint_vec_lbC_.resize(ncs+nch);
            constraint_vec_lbC_.segment(0,ncs)   = lds;
            constraint_vec_lbC_.segment(ncs,nch) = ldh;

            constraint_vec_ubC_.resize(ncs+nch);
            constraint_vec_ubC_.segment(0,ncs)   = uds;
            constraint_vec_ubC_.segment(ncs,nch) = udh;

        } else {
            constraint_mat_.resize(ncs+nch + nc_pre, dim_var_+ncs);
            constraint_mat_.setZero();
            constraint_mat_.block(0,0, ncs, dim_var_) = Cis;
            constraint_mat_.block(0,dim_var_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
            constraint_mat_.block(ncs,0, nch, dim_var_) = Cih;
            constraint_mat_.block(ncs+nch,0, nc_pre, dim_var_) = C_;

            constraint_vec_lbC_.resize(ncs+nch+nc_pre);
            constraint_vec_lbC_.segment(0,ncs)           = lds;
            constraint_vec_lbC_.segment(ncs,nch)         = ldh;
            constraint_vec_lbC_.segment(ncs+nch, nc_pre) = C_ld_;

            constraint_vec_ubC_.resize(ncs+nch+nc_pre);
            constraint_vec_ubC_.segment(0,ncs)          = uds;
            constraint_vec_ubC_.segment(ncs,nch)        = udh;
            constraint_vec_ubC_.segment(ncs+nch,nc_pre) = C_ud_;
        }

        // store new Constraints will be used next level
        C_.conservativeResize(nc_pre+nch+nt+ncs,dim_var_);
        C_.block(nc_pre,0, nch, dim_var_)        = Cih;
        C_.block(nc_pre+nch,0, nt, dim_var_)     = Ai;
        C_.block(nc_pre+nch+nt,0, ncs, dim_var_) = Cis;
        C_ld_.conservativeResize(nc_pre + nch + nt + ncs);
        C_ld_.segment(nc_pre, nch)        = ldh;
        C_ld_.segment(nc_pre+nch, nt)     = Eigen::VectorXd::Zero(nt);
        C_ld_.segment(nc_pre+nch+nt, ncs) = lds;
        C_ud_.conservativeResize(nc_pre + nch + nt + ncs);
        C_ud_.segment(nc_pre, nch)        = udh;
        C_ud_.segment(nc_pre+nch, nt)     = Eigen::VectorXd::Zero(nt);
        C_ud_.segment(nc_pre+nch+nt, ncs) = uds;

    } else if ( options_.hqp_algorithm_ == HQP_Nullspace ) {
        dim_u_ = NullBase_.cols();

        if (nv_i_ != (dim_u_ + ncs)){
            nv_i_ = dim_u_ + ncs;
            newQP = true;
        }
        // task mat
        task_mat_.resize(nt+ncs, nv_i_);
        task_mat_.setZero();
        task_mat_.block(0,0, nt, dim_u_) = wi.asDiagonal()*Ai*NullBase_;
        // task_mat_.block(nt, dim_u_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
        task_mat_.block(nt, dim_u_, ncs, ncs) = wis.asDiagonal();

        target_vec_.resize(nt+ncs);
        target_vec_.setZero();
        target_vec_.segment(0, nt) =  wi.asDiagonal()*( bi - Ai*xstar_ );

        int nc_pre = C_.rows();
        if (nc_i_ != (ncs + nch + nc_pre)) {
            nc_i_ = ncs + nch + nc_pre;
            newQP = true;
        }
        // constraint
        if (nc_pre == 0) {
            constraint_mat_.resize(ncs+nch, nv_i_);
            constraint_mat_.setZero();
            constraint_mat_.block(0,0, ncs, dim_u_) = Cis*NullBase_;
            constraint_mat_.block(0,dim_u_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
            constraint_mat_.block(ncs,0, nch, dim_u_) = Cih*NullBase_;

            constraint_vec_lbC_.resize(ncs+nch);
            constraint_vec_lbC_.segment(0,ncs)   = lds - Cis*xstar_;
            constraint_vec_lbC_.segment(ncs,nch) = ldh - Cih*xstar_;

            constraint_vec_ubC_.resize(ncs+nch);
            constraint_vec_ubC_.segment(0,ncs)   = uds - Cis*xstar_;
            constraint_vec_ubC_.segment(ncs,nch) = udh - Cih*xstar_;

        } else {
            constraint_mat_.resize(ncs+nch + nc_pre, nv_i_);
            constraint_mat_.setZero();
            constraint_mat_.block(0,0, ncs, dim_u_) = Cis*NullBase_;
            constraint_mat_.block(0,dim_u_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
            constraint_mat_.block(ncs,0, nch, dim_u_) = Cih*NullBase_;
            constraint_mat_.block(ncs+nch,0, nc_pre, dim_u_) = C_*NullBase_;

            constraint_vec_lbC_.resize(ncs+nch+nc_pre);
            constraint_vec_lbC_.segment(0,ncs)           = lds - Cis*xstar_;
            constraint_vec_lbC_.segment(ncs,nch)         = ldh - Cih*xstar_;
            constraint_vec_lbC_.segment(ncs+nch, nc_pre) = C_ld_ - C_*xstar_;

            constraint_vec_ubC_.resize(ncs+nch+nc_pre);
            constraint_vec_ubC_.segment(0,ncs)          = uds - Cis*xstar_;
            constraint_vec_ubC_.segment(ncs,nch)        = udh - Cih*xstar_;
            constraint_vec_ubC_.segment(ncs+nch,nc_pre) = C_ud_ - C_*xstar_;
        }

        // store new Constraints will be used next level
        C_.conservativeResize(nc_pre+nch+ncs,dim_var_);
        C_.block(nc_pre,0, nch, dim_var_)        = Cih;
        C_.block(nc_pre+nch,0, ncs, dim_var_) = Cis;
        C_ld_.conservativeResize(nc_pre + nch + ncs);
        C_ld_.segment(nc_pre, nch)        = ldh;
        C_ld_.segment(nc_pre+nch, ncs) = lds;
        C_ud_.conservativeResize(nc_pre + nch + ncs);
        C_ud_.segment(nc_pre, nch)        = udh;
        C_ud_.segment(nc_pre+nch, ncs) = uds;

    } else {
        assert(WBCKITS_ERROR && "wrong options parameter.");
    }
    
    // get QP
    hessian_mat_ = task_mat_.transpose()*task_mat_;
    regulation_ = options_.enable_regulation_;
    epsilon_ = options_.weight_regulation_;
    if (regulation_) {
        hessian_mat_.block(0,0,dim_u_,dim_u_) =
        hessian_mat_.block(0,0,dim_u_,dim_u_) + epsilon_*Eigen::MatrixXd::Identity(dim_u_,dim_u_);
    }
    
    gradient_vec_ = -task_mat_.transpose()*target_vec_;

    constraint_mat_t_ = constraint_mat_.transpose();

    // // Update Bound
    // bound_vec_lbx_ = lbx_*Eigen::VectorXd::Ones(nv_i_);
    // bound_vec_ubx_ = ubx_*Eigen::VectorXd::Ones(nv_i_);

    if (i==1) newQP = true;

    // create new QP if need
    if (newQP) {
        qp_var_ = nv_i_;
        qp_con_ = nc_i_;
        createNewQPSolver();
    }

    return true;

}

bool HQP::solveQP(int i)
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
            std::cout << "HQP solve failed at priority level "<<i<<".\n"
                      << "Error Message:" << qpOASES::MessageHandling::getErrorCodeMessage(status_code_solving_) << std::endl;
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
            std::cout << "HQP solve failed at priority level "<<i<<".\n"
                      << "Error Message:" << qpOASES::MessageHandling::getErrorCodeMessage(status_code_solving_) << std::endl;
            qp_->reset();
            init_done_ = false;
            
            return false;
        }
    }
    qp_computation_cost_ += cpu_time_;
    return true;
}

// bool HQP::createNewQPSolver()
// {
//     delete qp_;
//     qp_ = new qpOASES::SQProblem(nv_i_, nc_i_, qpOASES::HST_SEMIDEF);

//     if ( !setQPOptions() ) return false;

//     init_done_ = false;

//     return true;
// }

void HQP::showQPInfo(int i)
{
    std::cout << "hessian:\n" << hessian_mat_ << std::endl
              << "gradient_vec :" << gradient_vec_.transpose() << std::endl
              << "target mat :\n" << task_mat_ << std::endl
              << "target vec : " << target_vec_.transpose() << std::endl;
    std::cout << "---------------| priority: "<<i<<" |-----------------" << std::endl;
    std::cout << "pre C content:\n" << C_ << std::endl
              << "lC : " << C_ld_.transpose() << std::endl
              << "uC : " << C_ud_.transpose() << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "xstar : " << xstar_.transpose() << std::endl;

    std::cout << "------------------------------------------------------" << std::endl;

    std::cout << "constraint_mat_:\n" << constraint_mat_t_.transpose() << std::endl
              << "constraint_vec_lbC_: " << constraint_vec_lbC_.transpose() << std::endl
              << "constraint_vec_ubC_: " << constraint_vec_ubC_.transpose() << std::endl
              << "bound_vec_lbx_:" << bound_vec_lbx_.transpose() << std::endl
              << "bound_vec_ubx_:" << bound_vec_ubx_.transpose() << std::endl;
    if ( options_.hqp_algorithm_ = HQP_Nullspace ) {
        std::cout << "------------------------------------------------------" << std::endl;
        std::cout << "Nullbase matrix:" << std::endl
                  << NullBase_ << std::endl;
    }
    
    std::cout << "--------------------- END ----------------------------" << std::endl;
    
}

END_NAMESPACE_WBCKITS
