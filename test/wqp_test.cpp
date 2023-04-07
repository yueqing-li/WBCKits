#include <WBCKits.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>


int main()
{
    printf("wqp test.\n");

    WBCKits::WQP wqp_wbc;
    WBCKits::Options options;
    // presettings
    options.get_slack_variable_ = false;
    wqp_wbc.SetOptions(options);

    int var_dim = 7;

    Eigen::MatrixXd A1, A2, A3;
    Eigen::VectorXd b1, b2, b3;
    Eigen::VectorXd w2, w3;

    Eigen::MatrixXd C1, C2;
    Eigen::VectorXd lb1, lb2, wc1;
    Eigen::VectorXd ub1, ub2;

    A1 = Eigen::MatrixXd::Zero(3, var_dim);
    A1.block(0,0,3,3) = Eigen::MatrixXd::Identity(3,3);
    b1 = Eigen::VectorXd::Ones(3);

    A2 = Eigen::MatrixXd::Zero(3, var_dim);
    A2.block(0,3,3,3) = Eigen::MatrixXd::Identity(3,3);
    b2 = Eigen::VectorXd::Zero(3);
    w2 = 3.0*Eigen::VectorXd::Ones(3);

    A3 = Eigen::MatrixXd::Zero(2, var_dim);
    A3.block(0,5,2,2) = Eigen::MatrixXd::Identity(2,2);
    b3 = 5*Eigen::VectorXd::Ones(2);
    w3 = 0.5*Eigen::VectorXd::Ones(2);

    C1 = Eigen::MatrixXd::Identity(var_dim,var_dim);
    lb1 = -3.0*Eigen::VectorXd::Ones(var_dim);
    ub1 =  3.0*Eigen::VectorXd::Ones(var_dim);
    wc1 = 10*Eigen::VectorXd::Ones(var_dim);

    C2 = Eigen::MatrixXd::Zero(2, var_dim);
    C2.block(0,0,2,2) = Eigen::MatrixXd::Identity(2,2);
    C2.block(0,2,2,2) = Eigen::MatrixXd::Identity(2,2);
    lb2.resize(2);
    ub2.resize(2);
    lb2 << -3.0, -2.0;
    ub2 << 3.0, 2.0;

    // Add tasks
    wqp_wbc.AddTask(A1, b1);  // default weight w = I
    wqp_wbc.AddTask(A2, b2, w2);
    wqp_wbc.AddTask(A3, b3, w3);

    // Add constraints
    wqp_wbc.AddConstraint(C1,lb1,ub1);
    wqp_wbc.AddConstraint(C2,lb2, ub2);  // default weight w = I

    
    // results
    Eigen::VectorXd result;
    if (wqp_wbc.SolveWBC()) {
        wqp_wbc.GetResult(result);
    }
    std::cout<<"result: "<< result.transpose() << std::endl;

    // wqp reset
    wqp_wbc.Reset();
    // Add tasks
    wqp_wbc.AddTask(A1, b1);  // default weight w = I
    wqp_wbc.AddTask(A2, b2);
    wqp_wbc.AddTask(A3, b3);

    // Add constraints
    wqp_wbc.AddConstraint(C1,lb1,ub1);  // 
    wqp_wbc.AddConstraint(C2,lb2,ub2);  // default weight w = I

    
    // results
    if (wqp_wbc.SolveWBC()) {
        wqp_wbc.GetResult(result);
    }

    std::cout<<"result(same weight): "<< result.transpose() << std::endl;

    return 0;
}