#pragma once
#include <vector>
#include <Eigen/Dense>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(Eigen::MatrixXd particles) = 0;

protected:
    class System* m_system = nullptr;
};

