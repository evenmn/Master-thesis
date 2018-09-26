#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class WaveFunction
{
private:
    int m_M;
    int m_N;
    int m_D;
    int m_norbitals;
    int m_sampling;
    double m_sigma_sqrd;
    double m_omega_sqrd;
public:
    WaveFunction() {}
    int setTrialWF(int N, int M, int D, int norbitals, int sampling, double sigma_sqrd, double omega);
    double Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v);
};

void matrix(const VectorXd &Xa, int O, int D, int P_half, MatrixXd &A);
