#include "energy.h"
//#include "basis.h"
#include "general_tools.h"
#include "eigen3/Eigen/Dense"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


int Energy::init(int N, int M, int D, int O, int sampling, double sigma_sqrd, double omega)
{
    m_N          = N;
    m_P          = M/D;
    m_Phalf      = M/(2*D);
    m_M          = M;
    m_D          = D;
    m_O          = O;
    m_sampling   = sampling;
    m_sigma_sqrd = sigma_sqrd;
    m_omega_sqrd = omega*omega;
}


void Energy::rij_element(const VectorXd &X, int i, int j, MatrixXd &Dist) {
    // Update element (i,j) in distance matrix

    double dist = 0;
    for(int d=0; d<m_D; d++) {
        double diff = X(m_D*i+d)-X(m_D*j+d);
        dist += diff*diff;
    }
    Dist(i,j) = 1/sqrt(dist);
}


void Energy::rij(const VectorXd &X, MatrixXd &Dist) {
    // Fill up the entire distance matrix

    for(int i=0; i<m_P; i++) {
        for(int j=0; j<i; j++) {
            Energy::rij_element(X, i, j, Dist);
        }
    }
}

void Energy::rij_cross(const VectorXd &X, int par, MatrixXd &Dist) {
    // Update Dist when particle par is changed

    // Update row
    for(int i=0; i<par; i++) {
        Energy::rij_element(X, par, i, Dist);
    }
    // Update column
    for(int i=par+1; i<m_P; i++) {
        Energy::rij_element(X, i, par, Dist);
    }
}


// === DERIVATIVE PART ===

double Energy::dH(double x, int n) {
    //Derivative of Hermite polynomial of n'th degree
    if(n == 0) {
        return 0;
    }
    else {
        return 2*n*H(x,n-1);
    }
}


/*
void derivative(const VectorXd &Xa, int O, int D, int k, MatrixXd &dA) {
    //Derivative of A matrix

    int M = Xa.size();                                // Number of free dimensions
    int P_half = int(M/D);                        // Number of particles

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    // Find relevant row
    int row = int(k/D);

    // Find indices of relevant row
    VectorXd a = VectorXd::Zero(D);
    int l = k%D;
    for(int i=0; i<D; i++) {
        a(i) = k-l+i;
    }

    // Find matrix
    for(int i=0; i<P_half; i++) {
        dA(row, i) = dH(Xa(k), order(i, l));
        for(int j=0; j<D; j++) {
            if(a(j) != k) {
                dA(row, i) *= H(Xa(a(j)), order(i, j));
            }
        }
    }
}


void derivative3(const VectorXd &Xa, int O, int D, int k, MatrixXd &dA) {
    //Derivative of A matrix

    int M = Xa.size();                                // Number of free dimensions
    int P_half = int(M/D);                        // Number of particles

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    // Find relevant row
    int row = int(k/D);

    // Find indices of relevant row
    VectorXd a = VectorXd::Zero(D);
    int l = k%D;
    for(int i=0; i<D; i++) {
        a(i) = k-l+i;
    }

    // Find matrix
    for(int i=0; i<P_half; i++) {
        dA(row, D*i+l) = dH(Xa(k), order(i, l));
        for(int j=0; j<D; j++) {
            if(a(j) != k) {
                dA(row, D*i+l) *= H(Xa(a(j)), order(i, j));
            }
        }
    }
}
*/

void Energy::dA_element(const VectorXd &Xa, int k, MatrixXd &dA) {
    //Update element of dA
    double h = 1;
}


void Energy::dA_row(const VectorXd &Xa, int k, MatrixXd &dA) {
    //Update row of dA

    MatrixXd order = MatrixXd::Zero(m_Phalf, m_D);
    list(m_O, m_D, order);

    // Find indices of relevant row
    VectorXd a = VectorXd::Zero(m_D);
    int l = k%m_D;
    for(int i=0; i<m_D; i++) {
        a(i) = k-l+i;
    }

    // Find matrix
    for(int i=0; i<m_Phalf; i++) {
        dA(k, i) = dH(Xa(k), order(i, l));
        for(int j=0; j<m_D; j++) {
            if(a(j) != k) {
                dA(k, i) *= H(Xa(a(j)), order(i, j));
            }
        }
    }
}


void Energy::dA_matrix(const VectorXd &Xa, MatrixXd &dA) {
    //Initialize the entire dA matrix

    for(int k=0; k<m_Phalf; k++) {
        Energy::dA_row(Xa, k, dA);
    }
}


double Energy::EL_calc(const VectorXd X, const VectorXd Xa, const VectorXd v, const MatrixXd W, const MatrixXd &Dist, const MatrixXd &A_up_inv, const MatrixXd &A_dn_inv, const MatrixXd &dA_up, const MatrixXd &dA_dn,\
                             int interaction, double &E_kin, double &E_ext, double &E_int) {
    /*Local energy calculations*/

    // Set parameters to zero
    E_kin = 0;          // Total kinetic energy
    E_ext = 0;          // Total external potetial energy
    E_int = 0;          // Total interaction energy

    // Declare Eigen vectors
    VectorXd e_n = VectorXd::Zero(m_N);
    VectorXd e_p = VectorXd::Zero(m_N);
    VectorXd diff = VectorXd::Zero(m_M);

    // Fill up vectors
    for(int i=0; i<m_N; i++) {
        e_n(i) = 1/(1 + exp(-v(i)));
        e_p(i) = 1/(1 + exp(v(i)));
    }

    for(int i=0; i<m_M; i++) {
        //diff(i) = energy(Xa, m_D, m_norbitals, i);


        for(int j=0; j<m_Phalf; j++) {
            if(i<m_M/2) {
                diff(i) += dA_up(i,j)*A_up_inv(j,int(i/m_D));
            }
            else {
                diff(i) += dA_dn(i-m_M/2,j)*A_dn_inv(j,int((i-m_M/2)/m_D));
            }
        }

    }

    // === ENERGY CALCULATION ===
    // Kinetic energy
    if(m_sampling==2) {
        E_kin += (W*e_n).transpose()*diff;
        E_kin += 0.25*(W.cwiseAbs2()*e_p.cwiseProduct(e_n)).sum();
        E_kin -= (1/(2*m_sigma_sqrd))*(Xa.transpose()*W)*e_n;
        E_kin += 0.5*((W.transpose()*W).cwiseProduct(e_n*e_n.transpose())).sum();

        E_kin -= 0.5*m_M * m_sigma_sqrd;
        E_kin += 0.25*Xa.transpose() * Xa;
        E_kin -= (double) (diff.transpose()*Xa);
        E_kin = -E_kin/(2 * m_sigma_sqrd * m_sigma_sqrd);
    }

    else {
        E_kin += 2*(W*e_n).transpose()*diff;
        E_kin += (W.cwiseAbs2()*e_p.cwiseProduct(e_n)).sum();
        E_kin -= (2/m_sigma_sqrd)*(Xa.transpose()*W)*e_n;
        E_kin += ((W.transpose()*W).cwiseProduct(e_n*e_n.transpose())).sum();

        E_kin -= m_M;
        E_kin += (double) (Xa.transpose() * Xa)/m_sigma_sqrd;
        E_kin -= 2*(double) (diff.transpose()*Xa);
        E_kin = -E_kin/(2 * m_sigma_sqrd);
    }

    // Interaction energy
    if(interaction) E_int = Dist.sum();

    // Harmonic oscillator potential
    E_ext = (double) (X.transpose() * X) * m_omega_sqrd/ 2;

    return E_kin + E_ext + E_int;
}



