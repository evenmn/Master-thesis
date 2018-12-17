#pragma once
#include <string>
#include "eigen3/Eigen/Dense"

int factorial(int n);
double binomial(int n, int p);
int orbitals();
std::string generate_filename(std::string name, std::string extension);
double H(double x, int n);
