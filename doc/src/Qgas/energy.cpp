#include "energy.h"
//#include "basis.h"
#include "general_tools.h"
#include "eigen3/Eigen/Dense"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;


void rij_element(const VectorXd &X, int D, int i, int j, MatrixXd &Dist) {
    // Update element (i,j) in distance matrix

    double dist = 0;
    for(int d=0; d<D; d++) {
        double diff = X(D*i+d)-X(D*j+d);
        dist += diff*diff;
    }
    Dist(i,j) = 1/sqrt(dist);
}


void rij(const VectorXd &X, int D, MatrixXd &Dist) {
    // Fill up the entire distance matrix

    int P = X.size()/D;

    for(int i=0; i<P; i++) {
        for(int j=0; j<i; j++) {
            rij_element(X, D, i, j, Dist);
        }
    }
}

void rij_cross(const VectorXd &X, int D, int par, MatrixXd &Dist) {
    // Update Dist when particle par is changed

    int P = X.size()/D;

    // Update row
    for(int i=0; i<par; i++) {
        rij_element(X, D, par, i, Dist);
    }
    // Update column
    for(int i=par+1; i<P; i++) {
        rij_element(X, D, i, par, Dist);
    }
}

double dH(double x, int n) {
    //Derivative of Hermite polynomial of n'th degree
    if(n == 0) {
        return 0;
    }
    else {
        return 2*n*H(x,n-1);
    }
}


// === DERIVATIVE PART ===


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


void derivative4(const VectorXd &Xa, int O, int D, int k, MatrixXd &dA) {
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
        dA(k, i) = dH(Xa(k), order(i, l));
        for(int j=0; j<D; j++) {
            if(a(j) != k) {
                dA(k, i) *= H(Xa(a(j)), order(i, j));
            }
        }
    }
}


void derivative5(const VectorXd &Xa, int O, int D, MatrixXd &dA) {
    //Derivative of A matrix

    int M = Xa.size();                                // Number of free dimensions

    for(int k=0; k<M/2; k++) {
        derivative4(Xa, O, D, k, dA);
    }
}


int Energy::init(int N, int M, int D, int norbitals, int sampling, double sigma_sqrd, double omega)
{
    m_N          = N;
    m_M          = M;
    m_D          = D;
    m_norbitals  = norbitals;
    m_sampling   = sampling;
    m_sigma_sqrd = sigma_sqrd;
    m_omega_sqrd = omega*omega;
}


void Deter(const VectorXd &Xa, VectorXd &diff) {
    // Determinant dependent part



    //int n_orbitals = magic_numbers_inverse(20)+1;

    //int P = 20;
    //int D = 2;

    /*
    MatrixXd D_up = MatrixXd::Ones(int(3),int(3));
    MatrixXd D_dn = MatrixXd::Ones(int(3),int(3));

    matrix(Xa.head(P), n_orbitals, D, D_up);
    matrix(Xa.tail(P), n_orbitals, D, D_dn);

    VectorXd X_up = VectorXd::Zero(2*P);
    VectorXd X_dn = VectorXd::Zero(2*P);
    for(int i=0; i<P; i++) {
        X_up(i) = X_up(i+P) = Xa(i);
        X_dn(i) = X_dn(i+P) = Xa(i+P);
    }

    for(int i=0; i<P; i++) {
        if(i % 2==0){
            diff(i) = 4*(X_up(i+3) - X_up(i+5))/D_up.determinant();
            diff(i+6) = 4*(X_dn(i+3) - X_dn(i+5))/D_dn.determinant();
        }
        else {
            diff(i) = 4*(X_up(i+3) - X_up(i+1))/D_up.determinant();
            diff(i+6) = 4*(X_dn(i+3) - X_dn(i+1))/D_dn.determinant();
        }
    }
    */

    //for(int i=0; i<2*P; i++) {
    //    diff(i) = energy(Xa, D, i);
    //}
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

    int P = m_M/m_D;

    for(int i=0; i<m_M; i++) {
        //diff(i) = energy(Xa, m_D, m_norbitals, i);


        for(int j=0; j<P/2; j++) {
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



