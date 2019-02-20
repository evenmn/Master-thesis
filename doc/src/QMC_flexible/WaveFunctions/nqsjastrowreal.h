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
    void updateArrays(Eigen::VectorXd positions, int pRand);
    void resetArrays();
    void initializeArrays(Eigen::VectorXd positions);
    void updateParameters(Eigen::MatrixXd parameters);
    double evaluate();
    double evaluateSqrd();
    double computeFirstDerivative(Eigen::VectorXd positions, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

private:
    int m_elementNumber = 1;
    int m_numberOfHiddenNodes = 0;
    double m_sigmaSqrd = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
};
