#include <WBCKits.hpp>
#include <eigen3/Eigen/Dense>
#include <iostream>


int main()
{

    int var_dim = 7;

    Eigen::MatrixXd A1, A2, A3;
    Eigen::VectorXd b1, b2, b3;
    Eigen::VectorXd w2, w3;  // task weight is not used in GHC

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

    printf("GHC test.\n");

    WBCKits::GHC ghc_wbc;
    WBCKits::Options options;
    options.nwsr_max_ = 1000;
    options.enable_regulation_ = true;
    options.Print();
    ghc_wbc.SetOptions(options);

    // Add tasks
    ghc_wbc.AddTask(A1,b1);
    ghc_wbc.AddTask(A2, b2, w2);
    ghc_wbc.AddTask(A3, b3, w3);

    // Add constraints
    ghc_wbc.AddConstraint(C1,lb1,ub1, true);
    ghc_wbc.AddConstraint(C2,lb2, ub2, false);

    //------ set priority similar to hqp
    Eigen::VectorXi phi;
    Eigen::VectorXi chi;
    phi.resize(3);
    chi.resize(2);
    phi << 0, 1, 2;
    chi << 1, 2;
    
    ghc_wbc.SetPriority(phi, chi);

    
    // results
    Eigen::VectorXd result;
    if (ghc_wbc.SolveWBC()) {
        ghc_wbc.GetResult(result);
    }

    std::cout<<"result A : "<< result.transpose() << std::endl;




    // ---------  reset GHC --------

    ghc_wbc.Reset();
    options.enable_regulation_ = true;
    ghc_wbc.SetOptions(options);
    // Add tasks
    ghc_wbc.AddTask(A1, b1);
    ghc_wbc.AddTask(A2, b2, w2);
    ghc_wbc.AddTask(A3, b3, w3);

    // Add constraints
    ghc_wbc.AddConstraint(C1,lb1,ub1);
    ghc_wbc.AddConstraint(C2,lb2, ub2, false);


    // In GHC, if all tasks have the same priorities, 
    // it equals to solve each task respectively and add them togethere.
    phi << 1, 1, 1;
    chi << 1, 1;

    ghc_wbc.SetPriority(phi, chi);

    // results
    if (ghc_wbc.SolveWBC()) {
        ghc_wbc.GetResult(result);
    }

    std::cout<<"result B: "<< result.transpose() << std::endl;

    return 0;
}