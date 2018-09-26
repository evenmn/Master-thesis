#include "wavefunction.h"
#include "general_tools.h"
#include "eigen3/Eigen/Dense"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


// === DETERMINANT MATRIX PART ===


double A_elements(const VectorXd &Xa, int P_half, int D, int O, int i, int j) {
    // Updates an element in A-matrix

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    double element = 1;
    for(int k=0; k<D; k++) {
        element *= H(Xa(D*i+k), order(j,k));
    }
    return element;
}

void A_rows(const VectorXd &Xa, int P_half, int D, int O, int i, MatrixXd &A) {
    // Updates a row in A-matrix

    for(int j=0; j<P_half; j++) {
        A(i,j) = A_elements(Xa, P_half, D, O, i, j);
    }
}


void matrix(const VectorXd &Xa, int O, int D, int P_half, MatrixXd &A) {
    // Update the entire matrix

    for(int j=0; j<P_half; j++) {
        A_rows(Xa, P_half, D, O, j, A);
    }
}

double Jastrow_NQS(VectorXd v) {
    //Neural Quantum State Wavefunction (NQS-WF)

    int N = v.size();
    double prod = 1;
    for(int i=0; i < N; i++) {
        prod *= (1 + exp(v(i)));
    }

    return prod;
}


double Gauss_WF(VectorXd Xa, double sigma_sqrd) {
    //Gaussian WF

    return exp(-(double)(Xa.transpose() * Xa)/(2 * sigma_sqrd));
}


double Slater(int D, int O, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd) {
    // Setting up Slater determinant

    int M = Xa.size();                                // Number of free dimensions
    int P_half = int(M/(2*D));                        // Number of particles

    MatrixXd D_up = MatrixXd::Ones(P_half,P_half);
    MatrixXd D_dn = MatrixXd::Ones(P_half,P_half);

    matrix(Xa.head(M/2), O, D, P_half, D_up);
    matrix(Xa.tail(M/2), O, D, P_half, D_dn);


    return D_up.determinant()*D_dn.determinant()*Gauss_WF(Xa, sigma_sqrd)*Jastrow_NQS(v);
}


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


double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v)
{
    //Unnormalized wave function

    double Prob = Slater(m_D, m_norbitals, Xa, v, m_sigma_sqrd);
    return Prob * Prob;
}

