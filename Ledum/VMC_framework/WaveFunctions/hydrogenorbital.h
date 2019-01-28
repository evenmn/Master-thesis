#pragma once
#include "wavefunction.h"

class HydrogenOrbital : public WaveFunction {
public:
    HydrogenOrbital(class System* system, double beta);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(Eigen::MatrixXd particles, int k);
    double computeDoubleDerivative(Eigen::MatrixXd particles);
    double computeFirstEnergyDerivative(Eigen::MatrixXd particles);
    double computeDoubleEnergyDerivative(Eigen::MatrixXd particles);
};
