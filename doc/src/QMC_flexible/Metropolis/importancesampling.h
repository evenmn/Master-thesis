#pragma once
#include "metropolis.h"
#include <Eigen/Dense>
#include <vector>

class ImportanceSampling : public Metropolis {
public:
    ImportanceSampling(System* system);
    bool acceptMove();
    double QuantumForce(const Eigen::VectorXd positions, int i);
    double GreenFuncSum(const Eigen::VectorXd newPositions);

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
    std::vector<class WaveFunction*>    m_waveFunctionVector;
};
