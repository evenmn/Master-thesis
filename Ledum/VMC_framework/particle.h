#pragma once
#include <vector>
#include <Eigen/Dense>

class Particle {
public:
    Particle();
    void setPosition(const Eigen::MatrixXd &position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    Eigen::MatrixXd getPosition() { return m_position; }
    //std::vector<double> getPosition() { return m_position; }

private:
    int     m_numberOfDimensions = 0;
    Eigen::MatrixXd m_position = Eigen::MatrixXd();
    //std::vector<double> m_position = std::vector<double>();
};

