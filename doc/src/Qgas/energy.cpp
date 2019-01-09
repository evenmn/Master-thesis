#include "energy.h"
#include "basis.h"
#include "general_tools.h"
#include "eigen3/Eigen/Dense"
#include "wavefunction.h"
#include "common.h"

#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;

void Energy::rij_element(int i, int j) {
    // Update element (i,j) in Dist_inv_invance matrix

    double dist = 0;
    for(int d=0; d<D; d++) {
        double diff = X(D*i+d)-X(D*j+d);
        dist += diff*diff;
    }
    Dist_inv(i,j) = 1/sqrt(dist);
    Dist(i,j) = sqrt(dist);
}


void Energy::rij() {
    // Fill up the entire distance matrix

    for(int i=0; i<P; i++) {
        for(int j=0; j<i; j++) {
            Energy::rij_element(i, j);
        }
    }
}

void Energy::rij_cross(int par) {
    // Update Dist_inv when particle par is changed

    // Update row
    for(int i=0; i<par; i++) {
        Energy::rij_element(par, i);
    }
    // Update column
    for(int i=par+1; i<P; i++) {
        Energy::rij_element(i, par);
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


double Energy::EL_calc(double &E_kin, double &E_ext, double &E_int) {
    //Local energy calculations

    // Set parameters to zero
    E_kin = 0;          // Total kinetic energy
    E_ext = 0;          // Total external potetial energy
    E_int = 0;          // Total interaction energy

    // Declare Eigen vectors
    diff = VectorXd::Zero(M);

    // Fill up vectors
    for(int i=0; i<N; i++) {
        e(i) = 1/(1 + exp(-v(i)));
        e_p(i) = 1/(1 + exp(v(i)));
    }

    for(int i=0; i<M; i++) {
        for(int j=0; j<P_half; j++) {
            if(i<M_half) {
                diff(i) += dA_up(i,j)*A_up_inv(j,int(i/2));
            }
            else {
                diff(i) += dA_dn(i-M_half,j)*A_dn_inv(j,int((i-M_half)/2));
            }
        }
    }

    Slater Slat;
    Jastrow Jast;

    // === ENERGY CALCULATION ===
    // Kinetic energy
    if(sampling==2) {
        E_kin += (W*e).transpose()*diff;
        E_kin += 0.25*(W.cwiseAbs2()*e_p.cwiseProduct(e)).sum();
        E_kin -= (1/(2*sigma_sqrd))*(Xa.transpose()*W)*e;
        E_kin += 0.5*((W.transpose()*W).cwiseProduct(e*e.transpose())).sum();

        E_kin -= 0.5*M * sigma_sqrd;
        E_kin += 0.25*Xa.transpose() * Xa;
        E_kin -= double(diff.transpose()*Xa);
        E_kin = -E_kin/(2 * sigma_sqrd * sigma_sqrd);
    }

    else {
        int engcal = 1;
        if(engcal==0) {
            E_kin += 2*(W*e).transpose()*diff;
            E_kin += (W.cwiseAbs2()*e_p.cwiseProduct(e)).sum();                     //NQS Jastrow secder
            E_kin -= (2/sigma_sqrd)*(Xa.transpose()*W)*e;
            E_kin += ((W.transpose()*W).cwiseProduct(e*e.transpose())).sum();

            E_kin -= M;                                                             //Gauss ML secder
            E_kin += double(Xa.transpose() * Xa) * omega * omega/sigma_sqrd;        //Xa^2/sigma^2
            E_kin -= double(2*diff.transpose()*Xa);
            E_kin = -E_kin/(2 * sigma_sqrd);
        }
        else if(engcal==1) {
            E_kin += Jast.Jastrow_NQS(v, 0, 2);
            E_kin += Jast.PadeJastrow(0, 2);
            E_kin += Slat.Gauss_ML(Xa, 0, 2);
            E_kin += Slat.SlaterDet(Xa, 0, 2);
            //E_kin += Slat.Gauss_partly(Xa, 0, 2);

            for(int k=0; k<M; k++) {
                double p1 = Jast.Jastrow_NQS(v, k, 1);
                double p2 = Slat.Gauss_ML(Xa, k, 1);
                double p3 = Slat.SlaterDet(Xa, k, 1);
                double p4 = Jast.PadeJastrow(k, 1);
                //double p5 = Slat.Gauss_partly(Xa, k, 1);
                double sum = p1 + p2 + p3 + p4;
                E_kin += sum*sum;
            }
            E_kin = -E_kin/(2);
        }
    }

    // Interaction energy
    if(interaction) E_int = Dist_inv.sum();

    // Harmonic oscillator potential
    E_ext = double(X.transpose() * X) * omega*omega/ 2;

    return E_kin + E_ext + E_int;
}
