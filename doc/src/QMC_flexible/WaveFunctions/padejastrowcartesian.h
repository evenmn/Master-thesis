#pragma once
#include "wavefunction.h"

class PadeJastrowCartesian : public WaveFunction {
public:
    PadeJastrowCartesian(class System* system, int elementNumber);
    double evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix);
    double computeFirstDerivative(const Eigen::VectorXd particles, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

    double f(int i, int j);
    double g(int i, int j, int k, int l);
    double beta(int i, int j);
    double gamma();

private:
    int m_elementNumber = 1;
};
