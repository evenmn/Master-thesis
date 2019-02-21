#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class GradientDescentMomentum : public Optimization {
public:
    GradientDescentMomentum(System* system, double gamma);
    Eigen::VectorXd getImmediateGradients(WaveFunction* waveFunction);
    Eigen::MatrixXd getAllImmediateGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    double m_gamma = 0;
    Eigen::MatrixXd m_v;
};
