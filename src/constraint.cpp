#include "WBCKits/constraint.h"
#include <cassert>

BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/
Constraint::Constraint()
{
    assert(WBCKITS_ERROR && "no constraint define here ");
}

Constraint::Constraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
                       const Eigen::VectorXd &ubC, bool hard_constraint)
{
    constr_dim_ = C.rows();
    var_dim_ = C.cols();
    assert(constr_dim_==lbC.size());
    assert(constr_dim_==ubC.size());

    C_ = C;
    lbC_ = lbC;
    ubC_ = ubC;

    w_ = Eigen::VectorXd::Ones(constr_dim_);

    hard_constraint_ = hard_constraint;

}

Constraint::Constraint(const Eigen::MatrixXd &C, const Eigen::VectorXd &lbC,
                       const Eigen::VectorXd &ubC, bool hard_constraint,
                       const Eigen::VectorXd &w)
{
    constr_dim_ = C.rows();
    var_dim_ = C.cols();
    assert(constr_dim_==lbC.size());
    assert(constr_dim_==ubC.size());

    C_ = C;
    lbC_ = lbC;
    ubC_ = ubC;

    assert(constr_dim_==w.size());
    w_ = w;

    hard_constraint_ = hard_constraint;
}

END_NAMESPACE_WBCKITS
