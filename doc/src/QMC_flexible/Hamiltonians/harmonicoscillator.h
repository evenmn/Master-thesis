#pragma once
#include "hamiltonian.h"
#include <vector>
#include <Eigen/Dense>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system);
    double computeLocalEnergy();

private:
    double m_omega = 0;
};

