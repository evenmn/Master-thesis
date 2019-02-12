#pragma once
#include <Eigen/Dense>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy() = 0;
    virtual ~Hamiltonian() = 0;

protected:
    class System* m_system = nullptr;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_parameters;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_interaction = false;
};

