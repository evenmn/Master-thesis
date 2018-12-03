#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class WaveFunction {
public:
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
    double Psi_value(int D, int O, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd);
    double Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v);
};

class Slater: public WaveFunction {
public:
    Slater() {}
    double A_elements(const VectorXd &Xa, int P_half, int D, int O, int i, int j);
    void A_rows(const VectorXd &Xa, int P_half, int D, int O, int i, MatrixXd &A);
    void matrix(const VectorXd &Xa, int O, int D, int P_half, MatrixXd &A);
    double Gauss_WF(VectorXd Xa, double sigma_sqrd);
    double SlaterDet(int D, int O, const VectorXd &Xa);
};

class Jastrow: public WaveFunction {
public:
    Jastrow() {}
    double Jastrow_NQS(const VectorXd &v);
};

void matrix(const VectorXd &Xa, int O, int D, int P_half, MatrixXd &A);
void A_rows(const VectorXd &Xa, int P_half, int D, int O, int i, MatrixXd &A);
