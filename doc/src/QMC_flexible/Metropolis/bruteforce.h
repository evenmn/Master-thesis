#pragma once
#include "metropolis.h"
#include <Eigen/Dense>

class BruteForce : public Metropolis {
public:
    BruteForce(System* system);
    bool acceptMove();

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
};
