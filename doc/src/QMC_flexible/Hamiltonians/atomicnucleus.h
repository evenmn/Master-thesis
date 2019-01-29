#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class AtomicNucleus : public Hamiltonian {
public:
    AtomicNucleus(System* system, double omega, int numberOfParticles, int numberOfDimensions);
    double computeLocalEnergy(Eigen::MatrixXd particles);

private:
    double m_omega = 0;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
};
