#pragma once
#include "wavefunction.h"

class NQSJastrowReal : public WaveFunction {
public:
    NQSJastrowReal(class System* system, int elementNumber);
    Eigen::VectorXd b();
    Eigen::MatrixXd W();
    Eigen::VectorXd v(Eigen::VectorXd positions);
    Eigen::VectorXd f(Eigen::VectorXd positions);
    Eigen::VectorXd g(Eigen::VectorXd positions);
    double evaluate(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(const Eigen::VectorXd positions, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

private:
    int m_elementNumber = 1;
    int m_numberOfHiddenNodes = 0;
    double m_sigmaSqrd = 1;
};
