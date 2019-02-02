#pragma once
#include "hamiltonian.h"
#include <Eigen/Dense>

class AtomicNucleus : public Hamiltonian {
public:
    AtomicNucleus(System* system, int Z);
    double computeLocalEnergy();

private:
    int m_Z = 1;
};
