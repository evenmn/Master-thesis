#include "eigen3/Eigen/Dense"
#include "gradient_descent.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


int GradientDescent::init(int sampling, double sigma_sqrd) {
    m_sampling   = sampling;
    m_sigma_sqrd = sigma_sqrd;
}

void GradientDescent::Gradient_a(const VectorXd &Xa, VectorXd &da) {

    if(m_sampling==2) {
        da = 0.5*Xa/m_sigma_sqrd;
    }
    else{
        da = Xa/m_sigma_sqrd;
    }
}

void GradientDescent::Gradient_b(const VectorXd &e, VectorXd &db) {

    if(m_sampling==2) {
        db = 0.5*e;
    }
    else{
        db = e;
    }
}

void GradientDescent::Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW) {

    if(m_sampling==2) {
        dW = 0.5*X*e.transpose()/m_sigma_sqrd;
    }
    else{
        dW = X*e.transpose()/m_sigma_sqrd;
    }
}
