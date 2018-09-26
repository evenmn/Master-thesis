#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Energy
{
private:
    int m_M;
    int m_N;
    int m_D;
    int m_O;
    int m_P;
    int m_Phalf;
    int m_sampling;
    double m_sigma_sqrd;
    double m_omega_sqrd;
public:
    Energy() {}
    int init(int N, int M, int D, int O, int sampling, double sigma_sqrd, double omega);
    double EL_calc(const VectorXd X, const VectorXd Xa, const VectorXd v, const MatrixXd W, const MatrixXd &Dist, const MatrixXd &A_up_inv, const MatrixXd &A_dn_inv, const MatrixXd &dA_up, const MatrixXd &dA_dn, int interaction, double &E_kin, double &E_ext, double &E_int);
    void dA_row(const VectorXd &Xa, int k, MatrixXd &dA);
    void dA_matrix(const VectorXd &Xa, MatrixXd &dA);
    void dA_element(const VectorXd &Xa, int k, MatrixXd &dA);
    double dH(double x, int n);
    void rij(const VectorXd &X, MatrixXd &Dist);
    void rij_cross(const VectorXd &X, int par, MatrixXd &Dist);
    void rij_element(const VectorXd &X, int i, int j, MatrixXd &Dist);
};
