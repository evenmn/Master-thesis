#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Energy
{
public:
    Energy() {}
    double EL_calc(const MatrixXd &Dist, const MatrixXd &A_up_inv, const MatrixXd &A_dn_inv, const MatrixXd &dA_up, const MatrixXd &dA_dn, double &E_kin, double &E_ext, double &E_int);
    void dA_row(const VectorXd &Xa, int k, MatrixXd &dA);
    void dA_matrix(const VectorXd &Xa, MatrixXd &dA);
    void dA_element(const VectorXd &Xa, int k, MatrixXd &dA);
    double dH(double x, int n);
    void rij(const VectorXd &X, MatrixXd &Dist);
    void rij_cross(const VectorXd &X, int par, MatrixXd &Dist);
    void rij_element(const VectorXd &X, int i, int j, MatrixXd &Dist);
};
