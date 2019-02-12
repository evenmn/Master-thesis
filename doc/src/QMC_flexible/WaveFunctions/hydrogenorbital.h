#pragma once
#include "wavefunction.h"

class HydrogenOrbital : public WaveFunction {
public:
    HydrogenOrbital(class System* system, double beta);
    double evaluate(Eigen::MatrixXd positions);
    double computeDerivative(Eigen::MatrixXd positions);
    double computeEnergyDerivative(Eigen::MatrixXd positions);
private:
    double m_beta = 0;
};
