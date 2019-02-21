#pragma once
#include "optimization.h"
#include <Eigen/Dense>

class BarzilaiBorwein : public Optimization {
public:
    BarzilaiBorwein(System* system);
    Eigen::VectorXd getImmediateGradients(WaveFunction* waveFunction);
    Eigen::MatrixXd getAllImmediateGradients();
    Eigen::MatrixXd updateParameters();
    Eigen::MatrixXd getEnergyGradient();

private:
    double m_omega = 0;
    double m_omega_sqrd = 0;
    Eigen::MatrixXd m_oldParameters;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_oldGradients;
    Eigen::MatrixXd m_gradients;
};
