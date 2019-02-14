#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class GradientDescent : public Optimization {
public:
    GradientDescent(System* system);
    Eigen::VectorXd getImmediateGradients(WaveFunction* waveFunction);
    Eigen::MatrixXd getAllImmediateGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
};
