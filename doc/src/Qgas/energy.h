#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Energy
{
public:
    Energy() {}
    double EL_calc(double &E_kin, double &E_ext, double &E_int);
    void dA_row(const VectorXd &Xa, int k, MatrixXd &dA);
    void dA_matrix(const VectorXd &Xa, MatrixXd &dA);
    void dA_element(const VectorXd &Xa, int k, MatrixXd &dA);
    double dH(double x, int n);
    void rij();
    void rij_cross(int par);
    void rij_element(int i, int j);
};
