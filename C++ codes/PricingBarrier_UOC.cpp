
#include "Pricing_headers.h"

using namespace std;
using namespace Eigen;

EulerScheme_BS::EulerScheme_BS(
int steps,
float maturity,
float start_value,
float interest_rate,
float volatility,
float barrier,
std::string scheme_type,
std::string barrier_type,
float barrier_param
):
   interest_rate_(interest_rate),
   start_value_(start_value),
   volatility_(volatility),
   steps_(steps),
   maturity_(maturity),
   barrier_(barrier),
   scheme_type_(scheme_type),
   barrier_type_(barrier_type),
   barrier_param_(barrier_param)
{

}

VectorXd EulerScheme_BS::Simulate_Scheme()
    {
        if(barrier_type_ == "S")
            {
                if(scheme_type_ == "D")
                    {
                        float step = maturity_ / float(steps_);
                        VectorXd discrete_path(steps_ + 1);

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path(0) = start_value_;

                        for (int j(1) ; j <= steps_ ; j++)
                            {
                                discrete_path(j) = discrete_path(j-1) + interest_rate_* discrete_path(j-1) * step + volatility_ * discrete_path(j-1) * sqrt(step) * d(gen);
                            }


                        return discrete_path;
                    }

                else
                    {
                        float step = maturity_ / float(steps_);
                        VectorXd discrete_path(steps_ + 1);
                        float weight_coeff(1);
                        float max_value(start_value_);

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path(0) = start_value_;

                        for (int j(1) ; j <= steps_ ; j++)
                            {

                                discrete_path(j) = discrete_path(j-1) * (1 + interest_rate_ * step + volatility_ * sqrt(step) * d(gen));

                                max_value = max(max_value, (float)discrete_path(j));

                                float A((discrete_path(j - 1) - barrier_) * (discrete_path(j) - barrier_));
                                float c;

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) / maturity_) * (A / (pow(discrete_path(j - 1) * volatility_, 2))));
                                    }

                                weight_coeff = weight_coeff * (1 - exp(-c));
                            }

                        VectorXd result(3);

                        result(0) = discrete_path(last);
                        result(1) = max_value;
                        result(2) = weight_coeff;

                        return result;

                    }
            }

        else

            {
                if(scheme_type_ == "D")
                    {
                        float step = maturity_ / float(steps_);
                        VectorXd discrete_path(steps_ + 1);

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path(0) = start_value_ / float(barrier_param_);

                        for (int j(1) ; j <= steps_ ; j++)
                            {
                                discrete_path(j) = discrete_path(j-1) * (1 + (interest_rate_ - (float(1) / (barrier_param_ + step * (j-1)))) * step + volatility_ * sqrt(step) * (d(gen))); ////
                            }

                        return discrete_path;
                    }
                else
                    {
                        float step = maturity_ / float(steps_);
                        VectorXd discrete_path(steps_ + 1);
                        float weight_coeff(1);
                        float max_value(start_value_ / barrier_param_);

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path(0) = start_value_ / barrier_param_;

                        for (int j(1) ; j <= steps_ ; j++)
                            {

                                discrete_path(j) = discrete_path(j-1) * (1 + (interest_rate_ - (float(1) / (barrier_param_ + step * (j-1)))) * step + volatility_ * sqrt(step) * d(gen));

                                max_value = max(max_value, (float)discrete_path(j));

                                float A((discrete_path(j - 1) - 1) * (discrete_path(j) - 1));
                                float c;

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) / maturity_) * (A / (pow(discrete_path(j - 1) * volatility_, 2))));
                                    }

                                weight_coeff = weight_coeff * (1 - exp(-c));
                            }

                        VectorXd result(3);

                        result(0) = discrete_path(last);
                        result(1) = max_value;
                        result(2) = weight_coeff;

                        return result;
                    }
            }


    }


//


float MonteCarloPricing::d_positive()
    {
        return (1.0 / float(volatility_ * sqrt(maturity_))) * log(start_value_ * exp((interest_rate_) * (maturity_)) / float(strike_)) + pow(2, -1) * volatility_ * sqrt(maturity_);
    }

float MonteCarloPricing::d_negative()
    {
        return (1.0 / float(volatility_ * sqrt(maturity_))) * log(start_value_ * exp((interest_rate_) * (maturity_)) / float(strike_)) - pow(2, -1) * volatility_ * sqrt(maturity_);
    }


//Black & Scholes formula for Vanilla Euro Call Option Pricing
float MonteCarloPricing::Call_Price_BS()
    {
        float d_pos(d_positive());
        float d_neg(d_negative());
        return start_value_ * 0.5 * erfc(-d_pos * M_SQRT1_2) - strike_ * exp(- interest_rate_ * (maturity_)) * 0.5 * erfc(-d_neg * M_SQRT1_2);
    }

//



MonteCarloPricing::MonteCarloPricing(int steps,
                                    float maturity,
                                    float start_value,
                                    float interest_rate,
                                    float volatility,
                                    int mc_sample_size,
                                    float CI_risk,
                                    float strike,
                                    float barrier,
                                    std::string scheme_type,
                                    std::string barrier_type,
                                    float barrier_param
                                    ):
                                        EulerScheme_BS(steps,
                                        maturity,
                                        start_value,
                                        interest_rate,
                                        volatility,
                                        barrier,
                                        scheme_type,
                                        barrier_type,
                                        barrier_param), ////
                                        mc_sample_size_(mc_sample_size),
                                        CI_risk_(CI_risk),
                                        strike_(strike)
{

}


VectorXd MonteCarloPricing::MC_estimation()
    {
        float MC(0);
        float MC2(0);

        if(barrier_type_ == "S")
            {
                if(scheme_type_ == "D")
                    {
                        for (int i(0) ; i < mc_sample_size_ ; i++)
                            {

                                VectorXd output;
                                output = Simulate_Scheme();

                                float max_value(output.maxCoeff());
                                float value_maturity(output(last));
                                if(max_value <= barrier_)
                                    {
                                        MC = MC + max(value_maturity - strike_, float(0));
                                        MC2 = MC2 + pow(max(value_maturity - strike_, float(0)), 2);
                                    }
                            }
                    }
                else
                    {
                       for (int i(0) ; i < mc_sample_size_ ; i++)
                            {

                                VectorXd output;
                                output = Simulate_Scheme();
                                float simulated_process_maturity(output(0));
                                float max_value(output(1));
                                float coeff(output(2));

                                float ind;
                                float spread(barrier_ - max_value);
                                if(spread < 0)
                                    {
                                        ind = 0;
                                    }
                                else
                                    {
                                        ind = 1;
                                    }

                                float value(max((float)simulated_process_maturity - strike_, float(0)) * ind * coeff);
                                MC = MC + value;
                                MC2 = MC2 + pow(value, 2);
                            }
                    }

            }

        else
            {
                if(scheme_type_ == "D")
                    {
                        for (int i(0) ; i < mc_sample_size_ ; i++)
                            {

                                VectorXd output;
                                output = Simulate_Scheme();

                                float max_value(output.maxCoeff());
                                float value_maturity(output(last));

                                if(max_value <= 1)
                                    {
                                        MC = MC + max((barrier_param_ + maturity_) * value_maturity - strike_, float(0));
                                        MC2 = MC2 + pow(max((barrier_param_ + maturity_) * value_maturity - strike_, float(0)), 2);
                                    }
                            }
                    }

                else
                    {
                        for (int i(0) ; i < mc_sample_size_ ; i++)
                            {

                                VectorXd output;
                                output = Simulate_Scheme();
                                float simulated_process_maturity(output(0));
                                float max_value(output(1));
                                float coeff(output(2));

                                float ind;
                                float spread(1 - max_value);
                                if(spread < 0)
                                    {
                                        ind = 0;
                                    }
                                else
                                    {
                                        ind = 1;
                                    }

                                float value(max(float((barrier_param_ + maturity_) * simulated_process_maturity) - strike_, float(0)) * ind * coeff);
                                MC = MC + value;
                                MC2 = MC2 + pow(value, 2);
                            }
                    }
            }

    MC = MC / float(mc_sample_size_);
    MC2 = MC2 / float(mc_sample_size_);

    double r(CI_risk_);
    double c(0.5 * erfc(-(1-r/float(2)) * M_SQRT1_2));

    float std_dev(sqrt(MC2 - pow(MC, 2)));
    VectorXd result(2);
    result(0) = Call_Price_BS() - exp(-interest_rate_) * MC;
    result(1) = 2 * c * std_dev * (1.0 / float(sqrt(mc_sample_size_)));

    return result;

    }



// R-R extrapolation

EulerScheme_BS_RR::EulerScheme_BS_RR(
int steps,
float maturity,
float start_value,
float interest_rate,
float volatility,
float barrier,
int RR_order,
std::string scheme_type
):
   interest_rate_(interest_rate),
   start_value_(start_value),
   volatility_(volatility),
   steps_(steps),
   maturity_(maturity),
   barrier_(barrier),
   RR_order_(RR_order),
   scheme_type_(scheme_type)
{

}



std::vector<vector<float>> EulerScheme_BS_RR::Simulate_Scheme_RR()
    {
        //float step(maturity_ / float(RR_order_ * steps_));

        if(scheme_type_ == "D")
            {
                if(RR_order_ == 2)
                    {

                        std::vector<float> discrete_path_1;
                        std::vector<float> discrete_path_2;

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path_1.push_back(start_value_);
                        discrete_path_2.push_back(start_value_);

                        VectorXd save_rv(RR_order_ * steps_);

                        for (int j(1) ; j <= RR_order_ * steps_ ; j++)
                            {
                                float step(maturity_ / float(RR_order_ * steps_));
                                float U(d(gen));
                                save_rv(j - 1) = U;

                                discrete_path_2.push_back(discrete_path_2[j-1] + interest_rate_* discrete_path_2[j-1] * step + volatility_ * discrete_path_2[j-1] * sqrt(step) * U);
                            }

                        for (int j(1) ; j <= steps_ ; j++)
                            {
                                float step(maturity_ / float((RR_order_-1) * steps_));
                                float U((save_rv(2 * j - 1) + save_rv(2 * j - 2)) * sqrt(2));

                                discrete_path_1.push_back(discrete_path_1[j-1] + interest_rate_* discrete_path_1[j-1] * step + volatility_ * discrete_path_1[j-1] * sqrt(step) * U);
                            }

                        std::vector<vector<float>> output;
                        output.push_back(discrete_path_1);
                        output.push_back(discrete_path_2);

                        return output;

                    }

                else
                    {

                        std::vector<float> discrete_path_1;
                        std::vector<float> discrete_path_2;
                        std::vector<float> discrete_path_3;

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path_1.push_back(start_value_);
                        discrete_path_2.push_back(start_value_);
                        discrete_path_3.push_back(start_value_);

                        for (int j(1) ; j <= steps_; j++)
                            {
                                float U_1(d(gen));
                                float U_2(d(gen));
                                float U_3(d(gen));
                                float U_4(d(gen));

                                float A_1(U_1);
                                float A_2((U_2 + U_3) / sqrt(2));
                                float A_3(U_4);

                                float step(maturity_ / float(RR_order_ * steps_));

                                discrete_path_3.push_back(discrete_path_3[j-1] + interest_rate_* discrete_path_3[j-1] * step + volatility_ * discrete_path_3[j-1] * sqrt(step) * A_1);
                                discrete_path_3.push_back(discrete_path_3[j-1] + interest_rate_* discrete_path_3[j-1] * step + volatility_ * discrete_path_3[j-1] * sqrt(step) * A_2);
                                discrete_path_3.push_back(discrete_path_3[j-1] + interest_rate_* discrete_path_3[j-1] * step + volatility_ * discrete_path_3[j-1] * sqrt(step) * A_3);

                                A_1 = (sqrt(2) * U_1 + U_2) / sqrt(3);
                                A_2 = (U_3 + sqrt(2) * U_4) / sqrt(3);

                                step = maturity_ / float((RR_order_-1) * steps_);

                                discrete_path_2.push_back(discrete_path_2[j-1] + interest_rate_* discrete_path_2[j-1] * step + volatility_ * discrete_path_2[j-1] * sqrt(step) * A_1);
                                discrete_path_2.push_back(discrete_path_2[j-1] + interest_rate_* discrete_path_2[j-1] * step + volatility_ * discrete_path_2[j-1] * sqrt(step) * A_2);


                                float A((A_1 + A_2) / sqrt(2));

                                step = maturity_ / float((RR_order_-2) * steps_);

                                discrete_path_1.push_back(discrete_path_1[j-1] + interest_rate_* discrete_path_1[j-1] * step + volatility_ * discrete_path_1[j-1] * sqrt(step) * A);
                            }

                        std::vector<vector<float>> output;
                        output.push_back(discrete_path_1);
                        output.push_back(discrete_path_2);
                        output.push_back(discrete_path_3);

                        return output;

                    }
            }

        else

            {
                if(RR_order_ == 2)
                    {
                        float weight_coeff_1(1);
                        float weight_coeff_2(1);
                        float max_value_1(start_value_);
                        float max_value_2(start_value_);


                        std::vector<float> discrete_path_1;
                        std::vector<float> discrete_path_2;

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path_1.push_back(start_value_);
                        discrete_path_2.push_back(start_value_);

                        VectorXd save_rv(RR_order_ * steps_);

                        for (int j(1) ; j <= RR_order_ * steps_ ; j++)
                            {

                                float U(d(gen));
                                save_rv(j - 1) = U;

                                float step(maturity_ / float(RR_order_ * steps_));

                                discrete_path_2.push_back(discrete_path_2[j - 1] * (1 + interest_rate_ * step + volatility_ * sqrt(step) * U));

                                max_value_2 = max(max_value_2, (float)discrete_path_2[j - 1]);

                                float A((discrete_path_2[j - 1] - barrier_) * (discrete_path_2[j] - barrier_));
                                float c;

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * RR_order_ / maturity_) * (A / (pow(discrete_path_2[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_2 = weight_coeff_2 * (1 - exp(-c));

                            }


                        for (int j(1) ; j <= steps_ ; j++)
                            {

                                float U((save_rv(2 * j - 1) + save_rv(2 * j - 2)) / sqrt(2));

                                float step(maturity_ / float((RR_order_-1) * steps_));

                                discrete_path_1.push_back(discrete_path_1[j - 1] * (1 + interest_rate_ * step + volatility_ * sqrt(step) * U));

                                max_value_1 = max(max_value_1, (float)discrete_path_1[j - 1]);

                                float A((discrete_path_1[j - 1] - barrier_) * (discrete_path_1[j] - barrier_));
                                float c;

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_ * (RR_order_ - 1)) / maturity_) * (A / (pow(discrete_path_1[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_1 = weight_coeff_1 * (1 - exp(-c));

                            }

                        std::vector<vector<float>> results;

                        std::vector<float> values_maturity;
                        values_maturity.push_back(discrete_path_1.back());
                        values_maturity.push_back(discrete_path_2.back());
                        results.push_back(values_maturity);

                        std::vector<float> values_max;
                        values_max.push_back(max_value_1);
                        values_max.push_back(max_value_2);
                        results.push_back(values_max);

                        std::vector<float> weight_coeffs;
                        weight_coeffs.push_back(weight_coeff_1);
                        weight_coeffs.push_back(weight_coeff_2);
                        results.push_back(weight_coeffs);

                        return results;

                    }

                else
                    {
                        float weight_coeff_1(1);
                        float weight_coeff_2(1);
                        float weight_coeff_3(1);
                        float max_value_1(start_value_);
                        float max_value_2(start_value_);
                        float max_value_3(start_value_);


                        std::vector<float> discrete_path_1;
                        std::vector<float> discrete_path_2;
                        std::vector<float> discrete_path_3;

                        std::random_device rd{};
                        std::mt19937 gen{rd()};

                        normal_distribution<> d{0,1};

                        discrete_path_1.push_back(start_value_);
                        discrete_path_2.push_back(start_value_);
                        discrete_path_3.push_back(start_value_);


                        for (int j(1) ; j <= steps_ ; j++)
                            {

                                float U_1(d(gen));
                                float U_2(d(gen));
                                float U_3(d(gen));
                                float U_4(d(gen));

                                float A_1(U_1);
                                float A_2((U_2 + U_3) / sqrt(2));
                                float A_3(U_4);

                                float step(maturity_ / float((RR_order_) * steps_));

                                discrete_path_3.push_back(discrete_path_3[j-1] * (1 + interest_rate_ * step + volatility_ * sqrt(step) * A_1));

                                max_value_3 = max(max_value_3, (float)discrete_path_3.back());

                                float A((discrete_path_3[j - 1] - barrier_) * (discrete_path_3[j] - barrier_));
                                float c;

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * RR_order_ / maturity_) * (A / (pow(discrete_path_3[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_3 = weight_coeff_3 * (1 - exp(-c));

                                ///

                                discrete_path_3.push_back(discrete_path_3[j - 1] * (1 + interest_rate_ * step + volatility_ * sqrt(step) * A_2));

                                max_value_3 = max(max_value_3, (float)discrete_path_3[j - 1]);

                                A = (discrete_path_3[j - 1] - barrier_) * (discrete_path_3[j] - barrier_);

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * RR_order_ / maturity_) * (A / (pow(discrete_path_3[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_3 = weight_coeff_3 * (1 - exp(-c));

                                ///

                                discrete_path_3.push_back(discrete_path_3[j - 1] * (1 + interest_rate_ * step + volatility_ * sqrt(step) * A_3));


                                max_value_3 = max(max_value_3, (float)discrete_path_3[j - 1]);

                                A = (discrete_path_3[j - 1] - barrier_) * (discrete_path_3[j] - barrier_);

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * RR_order_ / maturity_) * (A / (pow(discrete_path_3[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_3 = weight_coeff_3 * (1 - exp(-c));




                                A_1 = (sqrt(2) * U_1 + U_2) / sqrt(3);
                                A_2 = (U_3 + sqrt(2) * U_4) / sqrt(3);

                                step = (maturity_ / float((RR_order_ - 1) * steps_));


                                discrete_path_2.push_back(discrete_path_2[j-1] + interest_rate_* discrete_path_2[j-1] * step + volatility_ * discrete_path_2[j-1] * sqrt(step) * A_1);

                                max_value_2 = max(max_value_2, (float)discrete_path_2[j - 1]);

                                A = (discrete_path_2[j - 1] - barrier_) * (discrete_path_2[j] - barrier_);

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * (RR_order_ - 1) / maturity_) * (A / (pow(discrete_path_2[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_2 = weight_coeff_2 * (1 - exp(-c));

                                ///

                                discrete_path_2.push_back(discrete_path_2[j-1] + interest_rate_* discrete_path_2[j-1] * step + volatility_ * discrete_path_2[j-1] * sqrt(step) * A_2);

                                max_value_2 = max(max_value_2, (float)discrete_path_2[j - 1]);

                                A = (discrete_path_2[j - 1] - barrier_) * (discrete_path_2[j] - barrier_);

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * (RR_order_ - 1) / maturity_) * (A / (pow(discrete_path_2[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_2 = weight_coeff_2 * (1 - exp(-c));



                                A = (A_1 + A_2) / sqrt(2);

                                step = (maturity_ / float((RR_order_ - 2) * steps_));

                                discrete_path_1.push_back(discrete_path_1[j-1] + interest_rate_* discrete_path_1[j-1] * step + volatility_ * discrete_path_1[j-1] * sqrt(step) * A);

                                max_value_1 = max(max_value_1, (float)discrete_path_1.back());

                                A = (discrete_path_1[j - 1] - barrier_) * (discrete_path_1[j] - barrier_);

                                if(A < 0)
                                    {
                                        c = 0;
                                    }
                                else
                                    {
                                        c = ((2 * float(steps_) * (RR_order_ - 2) / maturity_) * (A / (pow(discrete_path_1[j - 1] * volatility_, 2))));
                                    }

                                weight_coeff_1 = weight_coeff_1 * (1 - exp(-c));



                            }

                        std::vector<vector<float>> results;

                        std::vector<float> values_maturity;
                        values_maturity.push_back(discrete_path_1.back());
                        values_maturity.push_back(discrete_path_2.back());
                        values_maturity.push_back(discrete_path_3.back());
                        results.push_back(values_maturity);

                        std::vector<float> values_max;
                        values_max.push_back(max_value_1);
                        values_max.push_back(max_value_2);
                        values_max.push_back(max_value_3);
                        results.push_back(values_max);

                        std::vector<float> weight_coeffs;
                        weight_coeffs.push_back(weight_coeff_1);
                        weight_coeffs.push_back(weight_coeff_2);
                        weight_coeffs.push_back(weight_coeff_3);
                        results.push_back(weight_coeffs);

                        return results;
                    }
            }

    }



RR_extrapolation::RR_extrapolation(int steps,
                                    float maturity,
                                    float start_value,
                                    float interest_rate,
                                    float volatility,
                                    int mc_sample_size,
                                    float CI_risk,
                                    float strike,
                                    float barrier,
                                    int RR_order,
                                    std::string scheme_type
                                    ):
                                        EulerScheme_BS_RR(steps,
                                        maturity,
                                        start_value,
                                        interest_rate,
                                        volatility,
                                        barrier,
                                        RR_order,
                                        scheme_type),
                                        mc_sample_size_(mc_sample_size),
                                        CI_risk_(CI_risk),
                                        strike_(strike)
{

}

//d+ & d- computation
float RR_extrapolation::d_positive()
    {
        return (1.0 / float(volatility_ * sqrt(maturity_))) * log(start_value_ * exp((interest_rate_) * (maturity_)) / float(strike_)) + pow(2, -1) * volatility_ * sqrt(maturity_);
    }

float RR_extrapolation::d_negative()
    {
        return (1.0 / float(volatility_ * sqrt(maturity_))) * log(start_value_ * exp((interest_rate_) * (maturity_)) / float(strike_)) - pow(2, -1) * volatility_ * sqrt(maturity_);
    }


//Black & Scholes formula for Vanilla Euro Call Option Pricing
float RR_extrapolation::Call_Price_BS()
    {
        float d_pos(d_positive());
        float d_neg(d_negative());
        return start_value_ * 0.5 * erfc(-d_pos * M_SQRT1_2) - strike_ * exp(- interest_rate_ * (maturity_)) * 0.5 * erfc(-d_neg * M_SQRT1_2);
    }

//



VectorXd RR_extrapolation::RR_estimation()
    {
        float MC(0);
        float MC2(0);
        VectorXd coeff(RR_order_);

        if(scheme_type_ == "D")
            {
                if(RR_order_ == 2)
                    {
                        coeff(0) = -(1 + sqrt(2));
                        coeff(1) = sqrt(2) * (1 + sqrt(2));
                    }
                else
                    {
                        coeff(0) = (sqrt(3) - sqrt(2)) / (2 * sqrt(2) - sqrt(3) - 1);
                        coeff(1) = -2 * ((sqrt(3) - 1) / (2 * sqrt(2) - sqrt(3) - 1));
                        coeff(2) = 3 * ((sqrt(2) - 1) / (2 * sqrt(2) - sqrt(3) -1));
                    }

                std::vector<std::vector<float>> output;

                        for (int i(0) ; i < mc_sample_size_ ; i++)
                            {
                                output = Simulate_Scheme_RR();
                                VectorXd max_values(RR_order_);
                                VectorXd values_maturity(RR_order_);
                                VectorXd payoffs(RR_order_);

                                for (int j(0) ; j < RR_order_ ; j++)
                                    {
                                        max_values(j) = *max_element(std::begin(output[j]), std::end(output[j]));
                                        values_maturity(j) = output[j].back();

                                        if(max_values(j) <= barrier_)
                                            {
                                                payoffs(j) = max(float(values_maturity(j)- strike_), float(0));
                                                //cout << payoffs(j) << endl;
                                            }
                                        else
                                            {
                                                payoffs(j) = 0;
                                            }
                                    }

                                MC = MC + coeff.transpose() * payoffs;
                                MC2 = MC2 + pow(coeff.transpose() * payoffs, 2);
                            }

            }

        else
            {
                if(RR_order_ == 2)
                    {
                        coeff(0) = -1;
                        coeff(1) = 2;
                    }
                else
                    {
                        coeff(0) = 1.0 / sqrt(2);
                        coeff(1) = -4;
                        coeff(2) = 9.0 / sqrt(2);
                    }

                std::vector<std::vector<float>> output;


                for (int i(0) ; i < mc_sample_size_ ; i++)
                    {
                        output = Simulate_Scheme_RR();

                        VectorXd payoffs(RR_order_);

                        for (int j(0) ; j < RR_order_ ; j++)
                            {
                                float max_value(output[1][j]);
                                float value_maturity(output[0][j]);
                                float coeff_weight(output[2][j]);


                                float ind;
                                float spread(barrier_ - max_value);
                                if(spread < 0)
                                    {
                                        ind = 0;
                                    }
                                else
                                    {
                                        ind = 1;
                                    }

                                payoffs[j] = max((float)value_maturity - strike_, float(0)) * ind * coeff_weight;

                            }

                        MC = MC + (coeff.transpose() * payoffs);
                        MC2 = MC2 + pow(coeff.transpose() * payoffs, 2);

                    }
            }



        MC = MC / float(mc_sample_size_);
        MC2 = MC2 / float(mc_sample_size_);
        double r(CI_risk_);
        double c(0.5 * erfc(-(1-r/float(2)) * M_SQRT1_2));

        float std_dev(sqrt(MC2 - pow(MC, 2)));
        VectorXd result(2);
        result(0) = Call_Price_BS() - exp(-interest_rate_) * MC;
        result(1) = 2 * c * std_dev * (1.0 / float(sqrt(mc_sample_size_)));

        return result;

    }
