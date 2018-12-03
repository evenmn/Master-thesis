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


// ADAM

void Optimization::ADAM_a(const double eta, int i, VectorXd &m, VectorXd &v, const double b1, const double b2, const double EL_avg, const double MC, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a) {
    // Calculating Gradient Descent a
    VectorXd DA = VectorXd::Zero(m_M);
    Total_Gradient_a(EL_avg, MC, daE_tot, da_tot, DA);
    m = m*b1 + (1-b1)*DA;
    v = v*b2 + (1-b2)*DA.cwiseAbs2();
    VectorXd m_ = m/(1-pow(b1, (i+1)));
    VectorXd v_ = v/(1-pow(b2, (i+1)));
    for(int j=0; j<m_M; j++) {
        opt_a(i) = m_(i)/sqrt(v_(i));
    }
}

void Optimization::ADAM_b(const double eta, int i, VectorXd &m, VectorXd &v, const double b1, const double b2, const double EL_avg, const double MC, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b) {
    // Calculating Gradient Descent a
    VectorXd DB = VectorXd::Zero(m_N);
    Total_Gradient_b(EL_avg, MC, dbE_tot, db_tot, DB);
    m = m*b1 + (1-b1)*DB;
    v = v*b2 + (1-b2)*DB.cwiseAbs2();
    VectorXd m_ = m/(1-pow(b1, (i+1)));
    VectorXd v_ = v/(1-pow(b2, (i+1)));
    for(int j=0; j<m_M; j++) {
        opt_b(i) = m_(i)/sqrt(v_(i));
    }
}

void Optimization::ADAM_W(const double eta, int i, MatrixXd &m, MatrixXd &v, const double b1, const double b2, const double EL_avg, const double MC, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W) {
    // Calculating Gradient Descent a
    MatrixXd DW = MatrixXd::Zero(m_M, m_N);
    Total_Gradient_W(EL_avg, MC, dWE_tot, dW_tot, DW);
    m = m*b1 + (1-b1)*DW;
    v = v*b2 + (1-b2)*DW.cwiseAbs2();
    MatrixXd m_ = m/(1-pow(b1, (i+1)));
    MatrixXd v_ = v/(1-pow(b2, (i+1)));
    for(int j=0; j<m_M; j++) {
        for(int k=0; k<m_N; k++) {
            opt_W(j,k) = m_(j,k)/sqrt(v_(j,k));
        }
    }
}
