#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, int elementNumber);
    double evaluate(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(int k);
    double computeSecondDerivative();
    void computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k);
    void computeSecondEnergyDerivative(Eigen::VectorXd &gradients);
private:
    int     m_elementNumber = 0;
    double  m_omega         = 0;
};
