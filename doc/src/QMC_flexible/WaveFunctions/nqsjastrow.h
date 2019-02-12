#pragma once
#include "wavefunction.h"

class NQSJastrow : public WaveFunction {
public:
    NQSJastrow(class System* system, int elementNumber);
    Eigen::VectorXd b();
    Eigen::MatrixXd W();
    Eigen::VectorXd v(Eigen::VectorXd particles);
    Eigen::VectorXd f(Eigen::VectorXd particles);
    Eigen::VectorXd g(Eigen::VectorXd particles);
    double evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(const Eigen::VectorXd particles, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

private:
    int m_elementNumber = 1;
    int m_numberOfHiddenNodes = 0;
    double m_sigmaSqrd = 1;
};
