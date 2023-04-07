#include "WBCKits/rhp.h"
#include <iostream>

BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/

RHP::RHP()
{
    options_.enable_regulation_ = true;
    return;
}

RHP::~RHP()
{
    ClearTasks();
    ClearConstraints();
    delete qp_;
    qp_ = NULL;
}

void RHP::SetPriority(const Eigen::MatrixXd &phi, const Eigen::VectorXi &chi)
{

    if (task_list_.size() != phi.cols()) {
        assert(WBCKITS_ERROR && "tasks num doesn't match list's size.");
    }
    if (constraint_list_.size() != chi.size()) {
        assert(WBCKITS_ERROR && "constraints num doesn't match list's size.");
    }

    phi_ = phi;
    chi_ = chi;
}

void RHP::SetPriority(const Eigen::VectorXi &phi, const Eigen::VectorXi &chi)
{
    if (task_list_.size() != phi.size()) {
        assert(WBCKITS_ERROR && "tasks num doesn't match list's size.");
    }
    if (constraint_list_.size() != chi.size()) {
        assert(WBCKITS_ERROR && "constraints num doesn't match list's size.");
    }

    Eigen::MatrixXd phi_mat;
    int nt = phi.size();
    int priority_num = phi.maxCoeff();
    phi_mat.resize(priority_num, nt);
    phi_mat.setOnes();
    
    for (int k=0; k<nt; k++) {
        for (int i=0; i<(phi(k)-1); i++) {
            phi_mat(i,k) = 0.;
        }
        if ( phi(k) == 0 ) {
            phi_mat.col(k) = Eigen::VectorXd::Zero(priority_num);
        }
    }

    phi_ = phi_mat;
    chi_ = chi;

}

bool RHP::SolveWBC()
{
    return solveRHP();
}

void RHP::GetResult(Eigen::VectorXd & optimal_result)
{
    optimal_result = xstar_;
}

/********************************
 *  P R I V A T E               *
 * ******************************/

void RHP::initRHP()
{
    xstar_ = Eigen::VectorXd::Zero(dim_var_);
    C_.resize(0, dim_var_);
    C_ld_.resize(0);
    C_ud_.resize(0);
    qp_computation_cost_ = 0.;
}

bool RHP::solveRHP()
{
    bool solve_flag = true;
    initRHP();
    int priority_level = phi_.rows();
    // std::cout << "priority_level:" <<priority_level <<std::endl;
    for (int i=1; i<=priority_level; i++) {
        getGeneralizedProjector(i-1);

        if (constructQP(i)) {
            if (solveQP(i)) {
                qpOASES::real_t *primal_ptr = new qpOASES::real_t[nv_i_];
                Eigen::VectorXd ustar, vstar, primal_opt;
                primal_opt.resize(nv_i_);
                qp_->getPrimalSolution(primal_ptr);
                for (int k=0; k<nv_i_; k++) {
                    primal_opt(k) = (double)(primal_ptr[k]);
                }
                delete primal_ptr;
                int ncs_i_new = nv_i_ - dim_var_;
                ustar = primal_opt.segment(0,dim_var_);
                vstar = primal_opt.segment(dim_var_, ncs_i_new);

                xstar_ = xstar_ + Projector_*ustar;

                // ShowQPInfo(i);

                if (ncs_i_new != 0) {
                    int nc = C_.rows();
                    C_ld_.segment(nc-ncs_i_new, ncs_i_new) 
                                = C_ld_.segment(nc-ncs_i_new, ncs_i_new) - vstar;
                    C_ud_.segment(nc-ncs_i_new, ncs_i_new)
                                = C_ud_.segment(nc-ncs_i_new, ncs_i_new) - vstar;
                }
            } else {
                // // solve failed this level, drop constraints in this level
                // int ncs_i_new = nv_i_ - dim_var_;
                // int nc = C_.rows();
                // C_.conservativeResize(nc-ncs_i_new, dim_var_);
                // C_ld_.conservativeResize(nc - ncs_i_new);
                // C_ud_.conservativeResize(nc - ncs_i_new);

                solve_flag = false;

            }
        }
    }
    return solve_flag;
}

bool RHP::selectAndSortTask(int i, Eigen::MatrixXd &Aisa, Eigen::VectorXd &bisa, Eigen::VectorXd &wisa, Eigen::VectorXd &Lambdai, double eps)
{
    // std::cout << "Select Task :"<<i<<std::endl;

    i = i-1;
    Eigen::VectorXi s;
    int nt = phi_.cols();
    s.resize(nt,1);

    // select tasks
    int ns = 0;  //num of selected tasks
    for (int k=0; k<nt; k++) {
        if (i==0) {
            if (abs(phi_(i,k)) > eps) { // phi_(i,k) neq 0
                s(ns) = k;
                ns++;
            }
        } else {
            if (abs(phi_(i,k)-phi_(i-1,k)) > eps) { // phi_(i,k) neq phi_(i-1,k)
                s(ns) = k;
                ns++;
            }
        }
    }

    s.conservativeResize(ns,1);


    // sort phi_i with descending order
    for (int j=0; j<ns; j++) {
        for (int k=0; k<ns-j-1; k++) {
            if (phi_(i,s(k)) < phi_(i,s(k+1))) {
                int tmp = s(k);
                s(k) = s(k+1);
                s(k+1) = tmp;
            }
        }
    }

    // construct A_i^{A,s} and Lambda_i
    Aisa.resize(dim_obj_,dim_var_);
    bisa.resize(dim_obj_);
    wisa.resize(dim_obj_);
    Lambdai.resize(dim_obj_);


    int startrows = 0;
    for (int j=0; j<ns; j++) {
        int dim = task_list_[s(j)]->TaskDim();
        Aisa.block(startrows,0,dim, dim_var_) = task_list_[s(j)]->GetMat();
        bisa.segment(startrows, dim) = task_list_[s(j)]->GetVec();
        //construct weight matrix of tasks in this level
        wisa.segment(startrows, dim) = task_list_[s(j)]->GetWet();
        Lambdai.segment(startrows, dim) = phi_(i,s(j))*Eigen::VectorXd::Ones(dim);
        startrows += dim;
    }
    Aisa.conservativeResize(startrows, dim_var_);
    bisa.conservativeResize(startrows,1);
    wisa.conservativeResize(startrows,1);
    Lambdai.conservativeResize(startrows,1);

    if (ns == 0) return false;

    return true;
}

void RHP::getOrthBasis(const Eigen::MatrixXd &Aisa,
                       Eigen::MatrixXd &B, Eigen::VectorXi &origin, int &r, double eps)
{
    Eigen::MatrixXd Js;
    Js = Aisa*Projector_;
    int n = Js.cols();
    int m = Js.rows();
    B.resize(m, n);
    origin.resize(m, 1);
    int i = 0;
    for (int k=0; k<m; k++) {
        if (i >= n) {  // full colum rank, no more calc need
            break;
        }
        B.row(i) = Js.row(k);
        for (int j=0; j<i; j++) {
            // orthogonalization
            B.row(i) = B.row(i) - (B.row(i)*(B.row(j).transpose()))*B.row(j);
        }
        float scalar = std::sqrt(B.row(i)*(B.row(i).transpose()));
        if (scalar > eps) {
            B.row(i) = B.row(i)/scalar;
            origin(i) = k;
            i++;
        }
    }
    r = i;
    B.conservativeResize(r, n);
    origin.conservativeResize(r, 1);
}

void RHP::getGeneralizedProjector(int i)
{
    // std::cout << "GetGP :"<<i<<std::endl;
    if (i==0) {
        Projector_ = Eigen::MatrixXd::Identity(dim_var_,dim_var_);
        return;
    }
    Eigen::MatrixXd Aisa;
    Eigen::VectorXd bisa;
    Eigen::VectorXd wisa;
    Eigen::VectorXd lambda;
    bool flag;
    flag = selectAndSortTask(i, Aisa, bisa,wisa, lambda);
    if (!flag) {
        return; // no change to projector in this level
    }

    Eigen::MatrixXd B;
    Eigen::VectorXi origin;
    int r;
    getOrthBasis(Aisa, B, origin, r);

    Eigen::VectorXd Lambda_r;
    Lambda_r.resize(r);
    for (int k=0; k<r; k++) Lambda_r(k) = lambda(origin(k));

    Eigen::MatrixXd P_i;
    P_i = Eigen::MatrixXd::Identity(dim_var_,dim_var_) - B.transpose()*Lambda_r.asDiagonal()*B;
    Projector_ = Projector_*P_i;
}

bool RHP::selectConstraint(int i, Eigen::MatrixXd &Cis, Eigen::VectorXd &lds, Eigen::VectorXd &uds,
            Eigen::VectorXd &wis, Eigen::MatrixXd &Cih, Eigen::VectorXd &ldh, Eigen::VectorXd &udh)
{
    // std::cout << "Select Constraint :"<<i<<std::endl;
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

bool RHP::constructQP(int i)
{
    bool newQP = false;

    Eigen::MatrixXd Ai;
    Eigen::VectorXd bi, wi, lambda;
    bool haveTask;
    haveTask = selectAndSortTask(i, Ai, bi, wi, lambda);

    Eigen::MatrixXd Cis, Cih;
    Eigen::VectorXd lds, uds, wis, ldh, udh;
    bool haveConstraint;
    
    haveConstraint  = selectConstraint(i, Cis, lds, uds, wis, Cih, ldh, udh);

    int ncs = 0;  // num of soft constraint
    int nch = 0;  // num of hard constraint
    int nt = 0;
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

    if (nv_i_ != (dim_var_ + ncs)){
        newQP = true;
    }
    // task mat
    nv_i_ = dim_var_ + ncs;
    task_mat_.resize(nt+ncs, nv_i_);
    task_mat_.setZero();
    task_mat_.block(0,0, nt, dim_var_) = wi.asDiagonal()*Ai*Projector_;
    // task_mat_.block(nt, dim_var_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
    task_mat_.block(nt, dim_var_, ncs, ncs) = wis.asDiagonal();

    target_vec_.resize(nt+ncs);
    target_vec_.setZero();
    target_vec_.segment(0, nt) =  wi.asDiagonal()*(lambda.asDiagonal()*(bi - Ai*xstar_));

    int nc_pre = C_.rows();
    if (nc_i_ != (ncs + nch + nc_pre)) {
        newQP = true;
    }
    // constraint
    nc_i_ = ncs + nch + nc_pre;
    if (nc_pre == 0) {
        constraint_mat_.resize(ncs+nch, dim_var_+ncs);
        constraint_mat_.setZero();
        constraint_mat_.block(0,0, ncs, dim_var_) = Cis*Projector_;
        constraint_mat_.block(0,dim_var_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
        constraint_mat_.block(ncs,0, nch, dim_var_) = Cih*Projector_;

        constraint_vec_lbC_.resize(ncs+nch);
        constraint_vec_lbC_.segment(0,ncs)   = lds - Cis*xstar_;
        constraint_vec_lbC_.segment(ncs,nch) = ldh - Cih*xstar_;

        constraint_vec_ubC_.resize(ncs+nch);
        constraint_vec_ubC_.segment(0,ncs)   = uds - Cis*xstar_;
        constraint_vec_ubC_.segment(ncs,nch) = udh - Cih*xstar_;

    } else {
        constraint_mat_.resize(ncs+nch + nc_pre, dim_var_+ncs);
        constraint_mat_.setZero();
        constraint_mat_.block(0,0, ncs, dim_var_) = Cis*Projector_;
        constraint_mat_.block(0,dim_var_, ncs, ncs) = Eigen::MatrixXd::Identity(ncs,ncs);
        constraint_mat_.block(ncs,0, nch, dim_var_) = Cih*Projector_;
        constraint_mat_.block(ncs+nch,0, nc_pre, dim_var_) = C_*Projector_;

        constraint_vec_lbC_.resize(ncs+nch+nc_pre);
        constraint_vec_lbC_.segment(0,ncs)           = lds - Cis*xstar_;
        constraint_vec_lbC_.segment(ncs,nch)         = ldh - Cih*xstar_;
        constraint_vec_lbC_.segment(ncs+nch, nc_pre) = C_ld_ - C_*xstar_;

        constraint_vec_ubC_.resize(ncs+nch+nc_pre);
        constraint_vec_ubC_.segment(0,ncs)          = uds - Cis*xstar_;
        constraint_vec_ubC_.segment(ncs,nch)        = udh - Cih*xstar_;
        constraint_vec_ubC_.segment(ncs+nch,nc_pre) = C_ud_ - C_*xstar_;
    }

    // store new Constraint if have
    C_.conservativeResize(nc_pre+nch+ncs,dim_var_);
    C_.block(nc_pre,0, nch, dim_var_)     = Cih;
    C_.block(nc_pre+nch,0, ncs, dim_var_) = Cis;
    C_ld_.conservativeResize(nc_pre + nch + ncs);
    C_ld_.segment(nc_pre, nch)     = ldh;
    C_ld_.segment(nc_pre+nch, ncs) = lds;
    C_ud_.conservativeResize(nc_pre + nch + ncs);
    C_ud_.segment(nc_pre, nch)     = udh;
    C_ud_.segment(nc_pre+nch, ncs) = uds;

    // get QP
    hessian_mat_ = task_mat_.transpose()*task_mat_;
    // force to regulate for safer solving
    regulation_ = options_.enable_regulation_;
    epsilon_ = options_.weight_regulation_;
    if (regulation_) {
        hessian_mat_.block(0,0,dim_var_,dim_var_) =
        hessian_mat_.block(0,0,dim_var_,dim_var_) + epsilon_*Eigen::MatrixXd::Identity(dim_var_,dim_var_);
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

bool RHP::solveQP(int i)
{
    // std::cout << "SolveQP :"<<i<<std::endl;

    if (!init_done_) {
        // nwsr_ = nwsr_des_;
        nwsr_ = options_.nwsr_max_;
        // cpu_time_ = cpu_time_des_;
        cpu_time_ = options_.cpu_time_max_;
        // showQPInfo(i);
        status_code_solving_ = qp_->init(hessian_mat_.data(), gradient_vec_.data(),
                                         constraint_mat_t_.data(),
                                         bound_vec_lbx_.data(), bound_vec_ubx_.data(),
                                         constraint_vec_lbC_.data(), constraint_vec_ubC_.data(),
                                         nwsr_, &cpu_time_);

        if (status_code_solving_ > 0) {
            std::cout << "RHP solve failed at priority level "<<i<<".\n"
                      << "Error Message:" << qpOASES::MessageHandling::getErrorCodeMessage(status_code_solving_) << std::endl;
            qp_->reset();
            
            init_done_ = false;
            return false;
        } else {
            init_done_ = true;
        }
    } else {
        // nwsr_ = nwsr_des_;
        nwsr_ = options_.nwsr_max_;
        // cpu_time_ = cpu_time_des_;
        cpu_time_ = options_.cpu_time_max_;
        status_code_solving_ = qp_->hotstart(hessian_mat_.data(), gradient_vec_.data(),
                                             constraint_mat_t_.data(),
                                             bound_vec_lbx_.data(), bound_vec_ubx_.data(),
                                             constraint_vec_lbC_.data(), constraint_vec_ubC_.data(),
                                             nwsr_, &cpu_time_);
        if (status_code_solving_ > 0) {
            std::cout << "RHP solve failed at priority level "<<i<<".\n"
                      << "Error Message:" << qpOASES::MessageHandling::getErrorCodeMessage(status_code_solving_) << std::endl;
            qp_->reset();
            init_done_ = false;
            
            return false;
        }
    }
    qp_computation_cost_ += cpu_time_;
    return true;
}

// bool RHP::createNewQPSolver()
// {
//     delete qp_;
//     qp_ = new qpOASES::SQProblem(nv_i_, nc_i_, qpOASES::HST_SEMIDEF);

//     if ( !setQPOptions() ) return false;

//     init_done_ = false;

//     return true;
// }

void RHP::showQPInfo(int i)
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
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "project matrix:" << std::endl
              << Projector_ << std::endl;
    std::cout << "--------------------- END ----------------------------" << std::endl;
    
}

END_NAMESPACE_WBCKITS