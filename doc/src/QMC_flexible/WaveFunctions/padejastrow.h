#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system, int elementNumber);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(int k);
    double computeSecondDerivative();
    double computeFirstEnergyDerivative();
    double computeSecondEnergyDerivative();

private:
    int m_elementNumber = 1;
};
