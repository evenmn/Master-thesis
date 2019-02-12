#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system, int elementNumber);
    double evaluate(Eigen::MatrixXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::MatrixXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(int k);
    double computeSecondDerivative();
    void computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k);
    void computeSecondEnergyDerivative(Eigen::VectorXd &gradients);

private:
    int m_elementNumber = 1;
};
