#pragma once
#include "optimization.h"

class GradientDescent : public Optimization {
public:
    GradientDescent(class System* system);
    void getEnergyGradient(double EL_avg, Eigen::MatrixXd grad_tot, Eigen::MatrixXd gradE_tot, Eigen::MatrixXd &gradients);
    /*double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeDoubleDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeDoubleEnergyDerivative(Eigen::MatrixXd particles); */
};
