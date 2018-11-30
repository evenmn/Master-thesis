#include "eigen3/Eigen/Dense"
#include "optimization.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


int Optimization::init(int sampling, double sigma_sqrd, int M, int N) {
    // Initialize Optimization
    m_sampling   = sampling;
    m_sigma_sqrd = sigma_sqrd;
    m_M = M;
    m_N = N;
}

void Optimization::Gradient_a(const VectorXd &Xa, VectorXd &da) {
    // Calculating d(Psi)/d(a)
    if(m_sampling==2) {
        da = 0.5*Xa/m_sigma_sqrd;
    }
    else{
        da = Xa/m_sigma_sqrd;
    }
}

void Optimization::Gradient_b(const VectorXd &e, VectorXd &db) {
    // Calculating d(Psi)/d(b)
    if(m_sampling==2) {
        db = 0.5*e;
    }
    else{
        db = e;
    }
}

void Optimization::Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW) {
    // Calculating d(Psi)/d(W)
    if(m_sampling==2) {
        dW = 0.5*X*e.transpose()/m_sigma_sqrd;
    }
    else{
        dW = X*e.transpose()/m_sigma_sqrd;
    }
}

void Optimization::Total_Gradient_a(const double EL_avg, const double MC, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &DA) {
    // Calculating d(E)/d(a)
    DA =  2*(daE_tot - EL_avg*da_tot)/MC;
}

void Optimization::Total_Gradient_b(const double EL_avg, const double MC, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &DB) {
    // Calculating d(E)/d(b)
    DB =  2*(dbE_tot - EL_avg*db_tot)/MC;
}

void Optimization::Total_Gradient_W(const double EL_avg, const double MC, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &DW) {
    // Calculating d(E)/d(W)
    DW =  2*(dWE_tot - EL_avg*dW_tot)/MC;
}


// GRADIENT DESCENT

void Optimization::GD_a(const double eta, const double EL_avg, const double MC, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a) {
    // Calculating Gradient Descent a
    VectorXd DA = VectorXd::Zero(m_M);
    Total_Gradient_a(EL_avg, MC, daE_tot, da_tot, DA);
    opt_a = eta*DA;
}

void Optimization::GD_b(const double eta, const double EL_avg, const double MC, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b) {
    // Calculating Gradient Descent b
    VectorXd DB = VectorXd::Zero(m_N);
    Total_Gradient_b(EL_avg, MC, dbE_tot, db_tot, DB);
    opt_b = eta*DB;
}

void Optimization::GD_W(const double eta, const double EL_avg, const double MC, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W) {
    // Calculating Gradient Descent W
    MatrixXd DW = MatrixXd::Zero(m_M, m_N);
    Total_Gradient_W(EL_avg, MC, dWE_tot, dW_tot, DW);
    opt_W = eta*DW;
}
