#include "wavefunction.h"
#include "general_tools.h"
#include "basis.h"
#include "eigen3/Eigen/Dense"
#include "common.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;

// === SLATER ===

double Slater::A_elements(const VectorXd &Xa, double f(double, int), int i, int j) {
    // Updates an element in A-matrix

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    double element = 1;
    for(int k=0; k<D; k++) {
        element *= f(Xa(D*i+k), order(j,k));
    }

    return element;
}

void Slater::A_rows(const VectorXd &Xa, double f(double, int), int i, MatrixXd &A) {
    // Updates a row in A-matrix

    for(int j=0; j<P_half; j++) {
        A(i,j) = A_elements(Xa, f, i, j);
    }
}

void Slater::matrix(const VectorXd &Xa, double f(double, int), MatrixXd &A) {
    // Update the entire matrix

    for(int j=0; j<P_half; j++) {
        A_rows(Xa, f, j, A);
    }
}

double Slater::Gauss_ML(const VectorXd &Xa) {
    // Biased Gaussian

    return exp(-(double)(Xa.transpose() * Xa)/(2 * sigma_sqrd));
}

double Slater::Gauss(const VectorXd &X, double alpha) {
    // Gaussian with a variational parameter

    return exp(-(double)(X.transpose() * X) * alpha);
}

double Slater::SlaterDet(const VectorXd &Xa, double f(double, int)) {
    // Setting up Slater determinant

    MatrixXd D_up = MatrixXd::Ones(P_half,P_half);
    MatrixXd D_dn = MatrixXd::Ones(P_half,P_half);

    matrix(Xa.head(M/2), f, D_up);
    matrix(Xa.tail(M/2), f, D_dn);

    return D_up.determinant()*D_dn.determinant();
}


// === JASTROW ===

double Jastrow::Jastrow_NQS(const VectorXd &v) {
    //Neural Quantum State Wavefunction (NQS-WF)

    double prod = 1;
    for(int i=0; i < N; i++) {
        prod *= (1 + exp(v(i)));
    }

    return prod;
}

/*
double Jastrow::PadeJastrow(const VectorXd &Xa) {
    //Pade-Jastrow factor
    double sum = 0;
    for(int i=0; i<m_M; i++) {
        for(int j=0; j<m_M; j++) {
            sum += 1;
        }
    }
}
*/

double WaveFunction::Psi_value(const VectorXd &Xa, const VectorXd &v) {
    // Setting up total wavefunction

    Slater Slat;
    Jastrow Jast;

    double slater = 0;
    double jastrow = 0;

    int system = 0;
    if(system == 0) {
        // Hermite functions and NQS
        slater = Slat.SlaterDet(Xa, H) * Slat.Gauss_ML(Xa);
        jastrow = Jast.Jastrow_NQS(v);
    }

    return slater * jastrow;
}

double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v) {
    //Unnormalized wave function

    double Prob = Psi_value(Xa, v);
    return Prob * Prob;
}

