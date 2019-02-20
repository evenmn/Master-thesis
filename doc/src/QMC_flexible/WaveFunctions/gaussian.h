#pragma once
#include "wavefunction.h"

class Gaussian : public WaveFunction {
public:
    Gaussian(class System* system, int elementNumber);
    void updateArrays(Eigen::VectorXd positions, int pRand);
    void resetArrays();
    void initializeArrays(Eigen::VectorXd positions);
    double evaluate(Eigen::VectorXd particles);
    double evaluateSqrd(Eigen::VectorXd particles);
    double computeFirstDerivative(const Eigen::VectorXd positions, int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(int k);
    Eigen::VectorXd computeSecondEnergyDerivative();
private:
    int     m_elementNumber = 0;
    double  m_omega         = 0;
    double  m_alpha         = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
};
