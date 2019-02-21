#pragma once
#include "wavefunction.h"

class NQSJastrow : public WaveFunction {
public:
    NQSJastrow(class System* system, int elementNumber);
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
    Eigen::MatrixXd m_W;
    Eigen::VectorXd m_b;
    Eigen::VectorXd m_oldV;
    Eigen::VectorXd m_v;
    Eigen::VectorXd m_oldN;
    Eigen::VectorXd m_n;
    Eigen::VectorXd m_oldP;
    Eigen::VectorXd m_p;
};
