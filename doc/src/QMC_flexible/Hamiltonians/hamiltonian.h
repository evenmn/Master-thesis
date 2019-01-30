#pragma once
#include <Eigen/Dense>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(Eigen::MatrixXd particles) = 0;

    virtual ~Hamiltonian() = 0;

protected:
    class System* m_system = nullptr;
};

