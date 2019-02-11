#pragma once
#include "wavefunction.h"

class NQSJastrow : public WaveFunction {
public:
    NQSJastrow(class System* system, int elementNumber);
    double evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(const Eigen::VectorXd particles, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

private:
    int m_elementNumber = 1;
};
