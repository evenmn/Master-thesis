#include "wavefunction.h"
#include "general_tools.h"
#include "basis.h"
#include "eigen3/Eigen/Dense"
#include "common.h"

#include <cstdarg>   // Variadic templates
#include <cmath>     // Math
#include <ctime>     // Time
#include <iostream>  // Function

using namespace std;
using namespace Eigen;

// === SLATER ===

double Slater::A_elements(const VectorXd &Xa, double f(double, int), int i, int j) {
    // Updates an element in A-matrix

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    double element = 1;
    for(int k=0; k<D; k++) {
        element *= f(sqrt(omega) * Xa(D*i+k), order(j,k));
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


double Slater::Gauss_ML(const VectorXd &Xa, int k, int type) {
    // Biased Gaussian
    if(type == 0) {
        return exp(-(double)(omega * Xa.transpose() * Xa)/(2 * sigma_sqrd));
    }
    else if(type == 1) {
        // First derivative
        return -double(omega * Xa(k))/sigma_sqrd;
    }
    else if(type == 2) {
        // Second derivative
        return -M/sigma_sqrd;
    }
}


double Slater::Gauss(const VectorXd &X, double alpha, int k, int type) {
    // Gaussian with a variational parameter

    if(type == 0) {
        return exp(-(double)(omega * X.transpose() * X) * alpha);
    }
    else if(type == 1) {
        // First derivative
        return double(omega * X(k))/alpha;
    }
    else if(type == 2) {
        // Second derivative
        return omega/alpha;
    }
}


double Slater::SlaterDet(const VectorXd &Xa, double f(double, int), int k, int type) {
    // Setting up Slater determinant

    if(type == 0) {
        MatrixXd D_up = MatrixXd::Ones(P_half,P_half);
        MatrixXd D_dn = MatrixXd::Ones(P_half,P_half);

        matrix(Xa.head(M_half), f, D_up);
        matrix(Xa.tail(M_half), f, D_dn);

        return D_up.determinant()*D_dn.determinant();
    }
    else if(type == 1) {
        // First derivative
        return diff(k);
    }
    else if(type == 2) {
        // Second derivative
        return -double(diff.transpose() * diff);
    }
}


// === JASTROW ===

double Jastrow::Jastrow_NQS(const VectorXd &v, int k, int type) {
    //Neural Quantum State Wavefunction (NQS-WF)

    if(type == 0) {
        double prod = 1;
        for(int i=0; i < N; i++) {
            prod *= (1 + exp(v(i)));
        }
        return prod;
    }
    else if(type == 1) {
        // First derivative
        return double(W.row(k) * e)/sigma_sqrd;
    }
    else if(type == 2) {
        // Second derivative
        return (W.cwiseAbs2()*e_p.cwiseProduct(e)).sum();
    }
}


double Jastrow::PadeJastrow() {
    //Pade-Jastrow factor
    double sum = 0;
    for(int i=0; i<M; i++) {
        for(int j=0; j<M; j++) {
            sum += Dist(i,j)/(1+Dist(i,j));
        }
    }
    return exp(sum);
}

double WTF(double f1, double f2, double f3) {
    // Returns total WF

    return f1 * f2 * f3;
}

double ENG(double f1(VectorXd, int, int), VectorXd &f11, double f2(VectorXd, int, int), VectorXd &f21, double f3(VectorXd, int, int), VectorXd &f31) {
    // Return kinetic energy

    double e1 = 0;
    for(int k=0; k<M; k++) {
        e1 += f1(f11, k, 1) + f2(f21, k, 1) + f3(f31, k, 1);
    }
    double e2 = f1(f11, 0, 2)*f1(f11, 0, 2) + f2(f21, 0, 2)*f2(f21, 0, 2) + f3(f31, 0, 2)*f3(f31, 0, 2);

    return e1*e1 + e2;
}

double WaveFunction::Psi_value(const VectorXd &Xa, const VectorXd &v) {
    // Setting up total wavefunction

    Slater Slat;
    Jastrow Jast;

    double slater = 0;
    double jastrow = 0;

    int system = 0;
    // Idea: Make function which takes an arbitrary number of functions as
    // arguments and return the total wf and kinetic energy
    double result;
    if(system == 0) {
        // Hermite functions and NQS
        //result = WTF(Slat.SlaterDet(Xa, H), Slat.Gauss_ML(Xa, 0, 0), Jast.Jastrow_NQS(v, 0, 0));


        slater = Slat.SlaterDet(Xa, H, 0, 0) * Slat.Gauss_ML(Xa, 0, 0);
        jastrow = Jast.Jastrow_NQS(v, 0, 0);

        //double e1 = Slat.Gauss_ML(Xa, 0, 2) * Slat.Gauss_ML(Xa, 0, 2) + Jast.Jastrow_NQS(v, 0, 2)*Jast.Jastrow_NQS(v, 0, 2);
        //double e2 = 0;
        //for(int i=0; i<M; i++) {
        //    e2 += Slat.Gauss_ML(Xa, i, 1) + Jast.Jastrow_NQS(v, i, 1);
        //}
        //double energy = e1 + e2*e2;

    }

    return slater * jastrow;
}

double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v) {
    //Unnormalized wave function

    double Prob = Psi_value(Xa, v);
    return Prob * Prob;
}

double WaveFunction::Psi_ratio() {
    //Ratio between new and old probability distribution

    return Psi_value_sqrd(X_newa, v_new)/Psi_value_sqrd(Xa, v);
}
