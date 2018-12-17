#include "energy.h"
#include "basis.h"
#include "general_tools.h"
#include "eigen3/Eigen/Dense"
#include "common.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


void Energy::rij_element(const VectorXd &X, int i, int j, MatrixXd &Dist) {
    // Update element (i,j) in distance matrix

    double dist = 0;
    for(int d=0; d<D; d++) {
        double diff = X(D*i+d)-X(D*j+d);
        dist += diff*diff;
    }
    Dist(i,j) = 1/sqrt(dist);
}


void Energy::rij(const VectorXd &X, MatrixXd &Dist) {
    // Fill up the entire distance matrix

    for(int i=0; i<P; i++) {
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
    for(int i=par+1; i<P; i++) {
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

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    // Find indices of relevant row
    VectorXd a = VectorXd::Zero(D);
    int l = k%D;
    for(int i=0; i<D; i++) {
        a(i) = k-l+i;
    }

    // Find matrix
    for(int i=0; i<P_half; i++) {
        dA(k, i) = dH(Xa(k), order(i, l));
        for(int j=0; j<D; j++) {
            if(a(j) != k) {
                dA(k, i) *= H(Xa(a(j)), order(i, j));
            }
        }
    }
}


void Energy::dA_matrix(const VectorXd &Xa, MatrixXd &dA) {
    //Initialize the entire dA matrix

    for(int k=0; k<P_half; k++) {
        Energy::dA_row(Xa, k, dA);
    }
}


double Energy::EL_calc(const MatrixXd &Dist, const MatrixXd &A_up_inv, const MatrixXd &A_dn_inv, const MatrixXd &dA_up, const MatrixXd &dA_dn, \
                             double &E_kin, double &E_ext, double &E_int) {
    /*Local energy calculations*/

    // Set parameters to zero
    E_kin = 0;          // Total kinetic energy
    E_ext = 0;          // Total external potetial energy
    E_int = 0;          // Total interaction energy

    // Declare Eigen vectors
    VectorXd e_n = VectorXd::Zero(N);
    VectorXd e_p = VectorXd::Zero(N);
    VectorXd diff = VectorXd::Zero(M);

    // Fill up vectors
    for(int i=0; i<N; i++) {
        e_n(i) = 1/(1 + exp(-v(i)));
        e_p(i) = 1/(1 + exp(v(i)));
    }

    for(int i=0; i<M; i++) {
        //diff(i) = energy(Xa, m_D, m_norbitals, i);


        for(int j=0; j<P_half; j++) {
            if(i<M/2) {
                diff(i) += dA_up(i,j)*A_up_inv(j,int(i/D));
            }
            else {
                diff(i) += dA_dn(i-M/2,j)*A_dn_inv(j,int((i-M/2)/D));
            }
        }

    }

    // === ENERGY CALCULATION ===
    // Kinetic energy
    if(sampling==2) {
        E_kin += (W*e_n).transpose()*diff;
        E_kin += 0.25*(W.cwiseAbs2()*e_p.cwiseProduct(e_n)).sum();
        E_kin -= (1/(2*sigma_sqrd))*(Xa.transpose()*W)*e_n;
        E_kin += 0.5*((W.transpose()*W).cwiseProduct(e_n*e_n.transpose())).sum();

        E_kin -= 0.5*M * sigma_sqrd;
        E_kin += 0.25*Xa.transpose() * Xa;
        E_kin -= (double) (diff.transpose()*Xa);
        E_kin = -E_kin/(2 * sigma_sqrd * sigma_sqrd);
    }

    else {
        E_kin += 2*(W*e_n).transpose()*diff;
        E_kin += (W.cwiseAbs2()*e_p.cwiseProduct(e_n)).sum();
        E_kin -= (2/sigma_sqrd)*(Xa.transpose()*W)*e_n;
        E_kin += ((W.transpose()*W).cwiseProduct(e_n*e_n.transpose())).sum();

        E_kin -= M;
        E_kin += (double) (Xa.transpose() * Xa)/sigma_sqrd;
        E_kin -= 2*(double) (diff.transpose()*Xa);
        E_kin = -E_kin/(2 * sigma_sqrd);
    }

    // Interaction energy
    if(interaction) E_int = Dist.sum();

    // Harmonic oscillator potential
    E_ext = (double) (X.transpose() * X) * omega*omega/ 2;

    return E_kin + E_ext + E_int;
}



