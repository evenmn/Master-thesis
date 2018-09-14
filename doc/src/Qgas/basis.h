#pragma once

#include "eigen3/Eigen/Dense"

using namespace Eigen;

double H(double x, int n);
void matrix(const VectorXd &Xa, int N, int D, MatrixXd &A);
double matrix_2(const VectorXd &Xa, MatrixXd &D_up, MatrixXd &D_dn);
double Slater(int D, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd);
double energy(const VectorXd &Xa, int N, int D, int k);
int factorial(int n);