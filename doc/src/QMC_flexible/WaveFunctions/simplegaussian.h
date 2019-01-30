#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeSecondDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeSecondEnergyDerivative(Eigen::MatrixXd particles);
private:
    double m_alpha = 0;
};
