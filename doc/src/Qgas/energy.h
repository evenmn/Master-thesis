#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Energy
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
    Energy() {}
    int init(int N, int M, int D, int norbitals, int sampling, double sigma_sqrd, double omega);
    double EL_calc(const VectorXd X, const VectorXd Xa, const VectorXd v, const MatrixXd W, const MatrixXd &Dist, const MatrixXd &A_up_inv, const MatrixXd &A_dn_inv, const MatrixXd &dA_up, const MatrixXd &dA_dn, int interaction, double &E_kin, double &E_ext, double &E_int);
};

void derivative3(const VectorXd &Xa, int O, int D, int k, MatrixXd &dA);
void derivative4(const VectorXd &Xa, int O, int D, int k, MatrixXd &dA);
void derivative5(const VectorXd &Xa, int O, int D, MatrixXd &dA);
void A_rows(const VectorXd &Xa, int P_half, int D, int O, int j, MatrixXd &A);
void rij(const VectorXd &X, int D, MatrixXd &Dist);
void rij_cross(const VectorXd &X, int D, int par, MatrixXd &Dist);
