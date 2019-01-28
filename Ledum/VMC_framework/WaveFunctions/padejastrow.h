#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system, double alpha);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeDoubleDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeDoubleEnergyDerivative(Eigen::MatrixXd particles);
};
