#pragma once
#include "hamiltonian.h"
#include <vector>
#include <Eigen/Dense>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega);
    //double computeLocalEnergy(std::vector<Particle*> particles);
    double computeLocalEnergy(Eigen::MatrixXd particles);

private:
    double m_omega = 0;
};

