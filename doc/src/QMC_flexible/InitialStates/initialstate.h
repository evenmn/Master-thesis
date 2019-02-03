#pragma once
#include <Eigen/Dense>

class InitialState {
public:
    InitialState(class System* system);
    virtual void setupInitialState() = 0;
    Eigen::MatrixXd getParticles() { return m_particles; }
    Eigen::MatrixXd getDistanceMatrix() { return m_distanceMatrix; }
    Eigen::VectorXd getRadialVector() { return m_radialVector; }

    virtual ~InitialState() = 0;

protected:
    class System* m_system = nullptr;
    int m_numberOfDimensions = 0;
    int m_numberOfParticles = 0;
    Eigen::MatrixXd m_particles;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::VectorXd m_radialVector;
};

