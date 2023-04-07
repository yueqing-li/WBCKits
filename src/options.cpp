#include "WBCKits/options.h"
#include <iostream>

BEGIN_NAMESPACE_WBCKITS

/********************************
 *  P U B L I C                 *
 * ******************************/
Options::Options()
{
    SetDefault();
}

Options::Options( const Options& rhs)
{
    copy( rhs );
}

Options::~Options()
{
}

Options& Options::operator=( const Options& rhs )
{
    if ( this != &rhs )
    {
        copy( rhs );
    }

    return *this;
}

bool Options::SetDefault()
{
    nwsr_max_                           = 100;
    cpu_time_max_                       = 1.0;

    enable_regulation_                  = false;
    weight_regulation_                  = 1e-6;

    get_slack_variable_                 = false;

    hqp_algorithm_                      = HQP_Nullspace;

    decompose_method_                   = Decompose_SVD;
    decompose_threshold_                = 1e-4;

    return true;
}

void Options::Print() const
{
    char hqp_algorithm_name[][9] = {"Origin", "Nullbase"};
    char decompose_name[][4]     = {"QR", "SVD", "COD"};
    char bool_value[][6]         = {"False", "True"};
    std::cout << "----------- WBCKits Options -------------" << std::endl
              << "<<-- QP solver settings -->>" << std::endl
              << "nwsr_max_             = " << nwsr_max_ << std::endl
              << "cpu_time_max_         = " << cpu_time_max_ << std::endl
              << "<<-- QP construcion -->>" << std::endl
              << "enable_regulation_    = " << bool_value[enable_regulation_] << std::endl
              << "weight_regulation_    = " << weight_regulation_ << std::endl
              << "<<--      HQP       -->>" << std::endl
              << "hqp_algorithm_        = " << hqp_algorithm_name[hqp_algorithm_] << std::endl
              << "decompose_method_     = " << decompose_name[decompose_method_] << std::endl
              << "decompose_threshold_  = " << decompose_threshold_ << std::endl
              << "<<--    Result      -->>" << std::endl
              << "get_slack_variable_   = " << bool_value[get_slack_variable_] << std::endl
              ;

}

/********************************
 *  P R O T E C T E D           *
 * ******************************/
bool Options::copy( const Options& rhs )
{
    nwsr_max_ = rhs.nwsr_max_;
    cpu_time_max_ = rhs.cpu_time_max_;

    enable_regulation_ = rhs.enable_regulation_;
    weight_regulation_ = rhs.weight_regulation_;

    hqp_algorithm_ = rhs.hqp_algorithm_;

    get_slack_variable_ = rhs.get_slack_variable_;

    decompose_method_ = rhs.decompose_method_;
    decompose_threshold_ = rhs.decompose_threshold_;

    return true;
}



END_NAMESPACE_WBCKITS