#pragma once
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class GradientDescent
{
private:
    int m_sampling;
    double m_sigma_sqrd;
public:
    GradientDescent() {}
    int init(int sampling, double sigma_sqrd);
    void Gradient_a(const VectorXd &a, VectorXd &da);
    void Gradient_b(const VectorXd &e, VectorXd &db);
    void Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW);
};
