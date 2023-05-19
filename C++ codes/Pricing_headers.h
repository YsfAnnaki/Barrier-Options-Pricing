#ifndef DEF_GBM
#define DEF_GBM

#include <iostream>
#include <string>
#include <vector>
//#include <armadillo>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/Cholesky>
#include <random>
#include <math.h>
//#include <boost/random/random_device.hpp>
//#include <boost/random/normal_distribution.hpp>

//using namespace arma;
using namespace Eigen;

//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(MatrixXf)


class EulerScheme_BS
{
    public:
        EulerScheme_BS();
        EulerScheme_BS(int steps, float maturity, float start_value, float interest_rate, float volatility, float barrier, std::string scheme_type, std::string barrier_type, float barrier_param);

        VectorXd Simulate_Scheme(); //Modif

    protected:
        int steps_;
        float interest_rate_;
        float volatility_;
        float maturity_;
        float start_value_;
        float barrier_;
        std::string scheme_type_;
        std::string barrier_type_;
        float barrier_param_;
};

class MonteCarloPricing: public EulerScheme_BS
{
    public:
        MonteCarloPricing();
        MonteCarloPricing(int steps, float maturity, float start_value, float interest_rate, float volatility, int mc_sample_size, float CI_risk, float strike, float barrier, std::string scheme_type, std::string barrier_type, float barrier_param = 0);

        VectorXd MC_estimation();
        float d_positive();
        float d_negative();
        float Call_Price_BS();

    protected:
        float mc_sample_size_;
        float CI_risk_;
        float strike_;

};


// R-R extrapolation

class EulerScheme_BS_RR
{
    public:
        EulerScheme_BS_RR();
        EulerScheme_BS_RR(int steps, float maturity, float start_value, float interest_rate, float volatility, float barrier, int RR_order, std::string scheme_type);

        std::vector<std::vector<float>> Simulate_Scheme_RR(); //Modif

    protected:
        int steps_;
        float interest_rate_;
        float volatility_;
        float maturity_;
        float start_value_;
        int RR_order_;
        float barrier_;
        std::string scheme_type_;

};

class RR_extrapolation: public EulerScheme_BS_RR
{
    public:
        RR_extrapolation();
        RR_extrapolation(int steps, float maturity, float start_value, float interest_rate, float volatility, int mc_sample_size, float CI_risk, float strike, float barrier, int RR_order, std::string scheme_type);

        VectorXd RR_estimation();
        float d_positive();
        float d_negative();
        float Call_Price_BS();

    protected:
        float mc_sample_size_;
        float CI_risk_;
        float strike_;
};








#endif

