#pragma once

#include "eigen3/Eigen/Dense"

using namespace Eigen;

double H(double x, int n);
double Slater(int D, int O, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd);
double energy(const VectorXd &Xa, int D, int O, int k);
void matrix(const VectorXd &Xa, int O, int D, int M_half, MatrixXd &A);
void derivative2(const VectorXd &Xa, int O, int D, MatrixXd &dA);
void derivative3(const VectorXd &Xa, int O, int D, int k, MatrixXd &dA);
void A_rows(const VectorXd &Xa, int M_half, int D, int O, int j, MatrixXd &A);
