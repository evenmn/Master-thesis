#pragma once
#include "optimization.h"

class GradientDescent : public Optimization {
public:
    GradientDescent(class System* system);
    double gradient(Eigen::MatrixXd particles);
    /*double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeDoubleDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeDoubleEnergyDerivative(Eigen::MatrixXd particles); */
};
