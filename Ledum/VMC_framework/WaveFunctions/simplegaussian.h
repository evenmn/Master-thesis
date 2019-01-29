#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(Eigen::MatrixXd particles);
    double computeDerivative(Eigen::MatrixXd particles);
    double computeEnergyDerivative(Eigen::MatrixXd particles);
};
