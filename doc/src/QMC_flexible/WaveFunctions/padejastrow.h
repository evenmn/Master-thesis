#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system, double alpha, double beta, Eigen::MatrixXd Gamma);
    double evaluate(Eigen::MatrixXd particles);
    double computeDerivative(Eigen::MatrixXd particles);
    double computeEnergyDerivative(Eigen::MatrixXd particles);

private:
    Eigen::MatrixXd m_Gamma;
};
