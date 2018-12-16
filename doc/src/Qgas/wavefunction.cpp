#include "wavefunction.h"
#include "general_tools.h"
#include "basis.h"
#include "eigen3/Eigen/Dense"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


int WaveFunction::setTrialWF(int N, int M, int D, int norbitals, int sampling, double sigma_sqrd, double omega)
{
    m_N          = N;
    m_M          = M;
    m_D          = D;
    m_norbitals  = norbitals;
    m_sampling   = sampling;
    m_sigma_sqrd = sigma_sqrd;
    m_omega_sqrd = omega*omega;
}

// === SLATER ===

double Slater::A_elements(const VectorXd &Xa, double f(double, int), int P_half, int D, int O, int i, int j) {
    // Updates an element in A-matrix

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    double element = 1;
    for(int k=0; k<D; k++) {
        element *= f(Xa(D*i+k), order(j,k));
    }

    return element;
}

void Slater::A_rows(const VectorXd &Xa, double f(double, int), int P_half, int D, int O, int i, MatrixXd &A) {
    // Updates a row in A-matrix

    for(int j=0; j<P_half; j++) {
        A(i,j) = A_elements(Xa, f, P_half, D, O, i, j);
    }
}

void Slater::matrix(const VectorXd &Xa, double f(double, int), int O, int D, int P_half, MatrixXd &A) {
    // Update the entire matrix

    for(int j=0; j<P_half; j++) {
        A_rows(Xa, f, P_half, D, O, j, A);
    }
}

double Slater::Gauss_ML(const VectorXd &Xa, double sigma_sqrd) {
    // Biased Gaussian

    return exp(-(double)(Xa.transpose() * Xa)/(2 * sigma_sqrd));
}

double Slater::Gauss(const VectorXd &X, double alpha) {
    // Gaussian with a variational parameter

    return exp(-(double)(X.transpose() * X) * alpha);
}

double Slater::SlaterDet(const VectorXd &Xa, double f(double, int), int D, int O) {
    // Setting up Slater determinant

    int M = Xa.size();                                // Number of free dimensions
    int P_half = int(M/(2*D));                        // Number of particles

    MatrixXd D_up = MatrixXd::Ones(P_half,P_half);
    MatrixXd D_dn = MatrixXd::Ones(P_half,P_half);

    matrix(Xa.head(M/2), f, O, D, P_half, D_up);
    matrix(Xa.tail(M/2), f, O, D, P_half, D_dn);

    return D_up.determinant()*D_dn.determinant();
}


// === JASTROW ===

double Jastrow::Jastrow_NQS(const VectorXd &v) {
    //Neural Quantum State Wavefunction (NQS-WF)

    int N = v.size();
    double prod = 1;
    for(int i=0; i < N; i++) {
        prod *= (1 + exp(v(i)));
    }

    return prod;
}


double WaveFunction::Psi_value(int D, int O, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd) {
    // Setting up total wavefunction

    Slater Slat;
    Jastrow Jast;

    double slater = 0;
    double jastrow = 0;

    int system = 0;
    if(system == 0) {
        // Hermite functions and NQS
        slater = Slat.SlaterDet(Xa, H, m_D, m_norbitals) * Slat.Gauss_ML(Xa, sigma_sqrd);
        jastrow = Jast.Jastrow_NQS(v);
    }

    return slater * jastrow;
}

double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v)
{
    //Unnormalized wave function

    double Prob = Psi_value(m_D, m_norbitals, Xa, v, m_sigma_sqrd);
    return Prob * Prob;
}

