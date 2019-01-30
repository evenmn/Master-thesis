#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system, double beta, Eigen::MatrixXd Gamma);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeSecondDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeSecondEnergyDerivative(Eigen::MatrixXd particles);

private:
    Eigen::MatrixXd m_Gamma;
    double m_beta = 0;
};
