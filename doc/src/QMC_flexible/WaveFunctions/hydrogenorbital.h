#pragma once
#include "wavefunction.h"

class HydrogenOrbital : public WaveFunction {
public:
    HydrogenOrbital(class System* system, double beta);
    double evaluate(Eigen::MatrixXd particles);
    double computeDerivative(Eigen::MatrixXd particles);
    double computeEnergyDerivative(Eigen::MatrixXd particles);
};
