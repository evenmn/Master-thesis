#pragma once
#include <Eigen/Dense>

class Metropolis {
public:
    Metropolis(class System* system);
    Eigen::VectorXd updatePositions()      { return m_positions; }
    Eigen::VectorXd updateRadialVector()   { return m_radialVector; }
    Eigen::MatrixXd updateDistanceMatrix() { return m_distanceMatrix; }
    virtual bool acceptMove() = 0;
    virtual ~Metropolis() = 0;

protected:
    class System* m_system = nullptr;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_distanceMatrix;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfFreeDimensions = 0;
    double m_stepLength = 0;
    double m_diff = 0.5;
};
