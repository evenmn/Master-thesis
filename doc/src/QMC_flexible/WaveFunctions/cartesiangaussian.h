#pragma once
#include "wavefunction.h"

class CartesianGaussian : public WaveFunction {
public:
    CartesianGaussian(class System* system, int elementNumber);
    double evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(int k);
    double computeSecondDerivative();
    void computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k);
    void computeSecondEnergyDerivative(Eigen::VectorXd &gradients);
private:
    int     m_elementNumber = 0;
    double  m_omega         = 0;
    double  m_alpha         = 0;
};
