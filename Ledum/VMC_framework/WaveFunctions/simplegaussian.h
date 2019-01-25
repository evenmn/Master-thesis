#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeDoubleDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeDoubleEnergyDerivative(Eigen::MatrixXd particles);
};
