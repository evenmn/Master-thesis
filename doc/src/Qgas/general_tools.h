#pragma once
#include <string>
#include "eigen3/Eigen/Dense"

int factorial(int n);
double binomial(int n, int p);
double orbitals(int P, int D);
std::string generate_filename(int sampling, int P, int D, int N, int MC, int interaction, double sigma, double omega, double eta, std::string name, std::string extension);
void list(int N, int D, Eigen::MatrixXd &order);
double H(double x, int n);
