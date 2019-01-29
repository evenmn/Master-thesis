#pragma once
#include <Eigen/Dense>

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    Eigen::MatrixXd getParticles() { return m_particles; }

protected:
    class System* m_system = nullptr;
    Eigen::MatrixXd m_particles;
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
};

