#pragma once
#include "hamiltonian.h"
#include <vector>
#include <Eigen/Dense>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, int numberOfParticles, int numberOfDimensions, bool interaction);
    double computeLocalEnergy(Eigen::MatrixXd particles);

private:
    double m_omega = 0;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_interaction = false;
};

