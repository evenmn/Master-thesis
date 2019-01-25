#pragma once
#include <vector>
#include <Eigen/Dense>

class Particle {
public:
    Particle();
    void setPosition(const Eigen::VectorXd &position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    Eigen::VectorXd getPosition() { return m_position; }

private:
    int     m_numberOfDimensions = 0;
    Eigen::VectorXd m_position = Eigen::VectorXd();
};

