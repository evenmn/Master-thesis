#include "eigen3/Eigen/Dense"
#include "optimization.h"
#include "common.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;

void Optimization::Gradient_a(const VectorXd &Xa, VectorXd &da) {
    // Calculating d(Psi)/d(a)
    if(sampling==2) {
        da = 0.5*Xa/sigma_sqrd;
    }
    else{
        da = Xa/sigma_sqrd;
    }
}

void Optimization::Gradient_b(const VectorXd &e, VectorXd &db) {
    // Calculating d(Psi)/d(b)
    if(sampling==2) {
        db = 0.5*e;
    }
    else{
        db = e;
    }
}

void Optimization::Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW) {
    // Calculating d(Psi)/d(W)
    if(sampling==2) {
        dW = 0.5*X*e.transpose()/sigma_sqrd;
    }
    else{
        dW = X*e.transpose()/sigma_sqrd;
    }
}

/*
double Optimization::Gradient_B() {
    // Calculating d(Psi)/d(B)
    double part1 = 0;
    for(int k=0; k<P; k++) {
        for(int j=0; j<k; j++) {
            double
            part1 += 0.5*A(k,j)
        }
    }
}
*/

// Total gradients

void Optimization::Total_Gradient_a(const double EL_avg, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &DA) {
    // Calculating d(E)/d(a)
    DA =  2*(daE_tot - EL_avg*da_tot)/MC;
}

void Optimization::Total_Gradient_b(const double EL_avg, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &DB) {
    // Calculating d(E)/d(b)
    DB =  2*(dbE_tot - EL_avg*db_tot)/MC;
}

void Optimization::Total_Gradient_W(const double EL_avg, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &DW) {
    // Calculating d(E)/d(W)
    DW =  2*(dWE_tot - EL_avg*dW_tot)/MC;
}


// GRADIENT DESCENT

void Optimization::GD_a(const double EL_avg, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a) {
    // Calculating Gradient Descent a
    VectorXd DA = VectorXd::Zero(M);
    Total_Gradient_a(EL_avg, daE_tot, da_tot, DA);
    opt_a = eta*DA;
}

void Optimization::GD_b(const double EL_avg, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b) {
    // Calculating Gradient Descent b
    VectorXd DB = VectorXd::Zero(N);
    Total_Gradient_b(EL_avg, dbE_tot, db_tot, DB);
    opt_b = eta*DB;
}

void Optimization::GD_W(const double EL_avg, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W) {
    // Calculating Gradient Descent W
    MatrixXd DW = MatrixXd::Zero(M, N);
    Total_Gradient_W(EL_avg, dWE_tot, dW_tot, DW);
    opt_W = eta*DW;
}

//void Optimization::GD_B(const )


// ADAM

void Optimization::ADAM_a(int i, VectorXd &m, VectorXd &v, const double b1, const double b2, const double EL_avg, const VectorXd &daE_tot, const VectorXd &da_tot, VectorXd &opt_a) {
    // Calculating Gradient Descent a
    VectorXd DA = VectorXd::Zero(M);
    Total_Gradient_a(EL_avg, daE_tot, da_tot, DA);
    m = m*b1 + (1-b1)*DA;
    v = v*b2 + (1-b2)*DA.cwiseAbs2();
    VectorXd m_ = m/(1-pow(b1, (i+1)));
    VectorXd v_ = v/(1-pow(b2, (i+1)));
    for(int j=0; j<M; j++) {
        opt_a(j) = eta*m_(j)/sqrt(v_(j));
    }
}

void Optimization::ADAM_b(int i, VectorXd &m, VectorXd &v, const double b1, const double b2, const double EL_avg, const VectorXd &dbE_tot, const VectorXd &db_tot, VectorXd &opt_b) {
    // Calculating Gradient Descent a
    VectorXd DB = VectorXd::Zero(N);
    Total_Gradient_b(EL_avg, dbE_tot, db_tot, DB);
    m = m*b1 + (1-b1)*DB;
    v = v*b2 + (1-b2)*DB.cwiseAbs2();
    VectorXd m_ = m/(1-pow(b1, (i+1)));
    VectorXd v_ = v/(1-pow(b2, (i+1)));
    for(int j=0; j<N; j++) {
        opt_b(j) = eta*m_(j)/sqrt(v_(j));
    }
}

void Optimization::ADAM_W(int i, MatrixXd &m, MatrixXd &v, const double b1, const double b2, const double EL_avg, const MatrixXd &dWE_tot, const MatrixXd &dW_tot, MatrixXd &opt_W) {
    // Calculating Gradient Descent a
    MatrixXd DW = MatrixXd::Zero(M, N);
    Total_Gradient_W(EL_avg, dWE_tot, dW_tot, DW);
    m = m*b1 + (1-b1)*DW;
    v = v*b2 + (1-b2)*DW.cwiseAbs2();
    MatrixXd m_ = m/(1-pow(b1, (i+1)));
    MatrixXd v_ = v/(1-pow(b2, (i+1)));
    for(int j=0; j<M; j++) {
        for(int k=0; k<N; k++) {
            opt_W(j,k) = eta*m_(j,k)/sqrt(v_(j,k));
        }
    }
}
