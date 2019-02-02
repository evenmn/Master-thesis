#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, int elementNumber);
    double evaluate(Eigen::MatrixXd particles);
    double computeFirstDerivative(int k);
    double computeSecondDerivative();
    double computeFirstEnergyDerivative();
    double computeSecondEnergyDerivative();
private:
    int m_elementNumber = 0;
};
