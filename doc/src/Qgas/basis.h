#pragma once

#include "eigen3/Eigen/Dense"

using namespace Eigen;

double H(double x, int n);
double Slater(int D, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd);
double energy(const VectorXd &Xa, int D, int k);
