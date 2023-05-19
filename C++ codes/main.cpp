#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <random>

#include "Pricing_headers.h"
#include <time.h>
#include <eigen3/Eigen/Eigenvalues>
//#include "matplotlibcpp.h"

using namespace std;
using namespace Eigen;

//namespace plt = matplotlibcpp;


int main()
{
    //Parameters:

    float interest_rate(0.01);
    float volatility(0.2);
    float start_value(100);
    float maturity(1);
    int steps(1000);
    int mc_sample_size(100000);
    float CI_risk(0.01);
    float barrier(120);
    float barrier_param(120);
    float strike(110);
    std::string barrier_type("S");
    std::string scheme_type("C");
    int RR_order(3);


    //Pricing of an UIC contract with static barrier: (Put into comment the following 4 lignes if not needed)

    MonteCarloPricing MC(steps, maturity, start_value, interest_rate, volatility, mc_sample_size, CI_risk, strike, barrier, scheme_type, barrier_type);
    VectorXd result;
    result = MC.MC_estimation();
    cout << result(0) << endl;  // Print the estimator


    //Pricing of an UIC contract with moving barrier: (Put into comment the following 4 lignes if not needed)

    MonteCarloPricing MC(steps, maturity, start_value, interest_rate, volatility, mc_sample_size, CI_risk, strike, barrier, scheme_type, barrier_type, barrier_param);
    VectorXd result;
    result = MC.MC_estimation();
    cout << result(0) << endl;  // Print the estimator


    //RR extrapolation: (Put into comment the following 4 lignes if not needed)

    RR_extrapolation MC(steps, maturity, start_value, interest_rate, volatility, mc_sample_size, CI_risk, strike, barrier, RR_order, scheme_type);

    VectorXd result;
    result = MC.RR_estimation();
    cout << result(0) << endl;  // Print the estimator

    return 0;

}


