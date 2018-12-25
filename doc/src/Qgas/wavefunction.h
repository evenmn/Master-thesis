#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class WaveFunction {
public:
    WaveFunction() {}
    double Psi_value(const VectorXd &Xa, const VectorXd &v);
    double Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v);
    double Psi_ratio();
};

class Slater: public WaveFunction {
public:
    Slater() {}
    double A_elements(const VectorXd &Xa, double f(double, int), int i, int j);
    void A_rows(const VectorXd &Xa, double f(double, int), int i, MatrixXd &A);
    void matrix(const VectorXd &Xa, double f(double, int), MatrixXd &A);
    double Gauss(const VectorXd &X, double alpha, int k, int type);
    //double Gauss_ML(const VectorXd &Xa);
    double Gauss_ML(const VectorXd &Xa, int k, int type);
    double SlaterDet(const VectorXd &Xa, double f(double, int), int k, int type);
};

class Jastrow: public WaveFunction {
public:
    Jastrow() {}
    double Jastrow_NQS(const VectorXd &v, int k, int type);
    double PadeJastrow();
};

void matrix(const VectorXd &Xa, int O, int D, int P_half, MatrixXd &A);
void A_rows(const VectorXd &Xa, int P_half, int D, int O, int i, MatrixXd &A);
