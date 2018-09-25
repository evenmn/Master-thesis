#include <iostream>
#include <ctime>
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "general_tools.h"
#include "basis.h"

using namespace Eigen;
using namespace std;


// === BASIS PART ===

double H(double x, int n) {
    //Hermite polynomial of n'th degree

    if(n == 0) {
        return 1;
    }
    else if(n == 1) {
        return 2*x;
    }
    else {
        return 2*x*H(x,n-1)-2*(n-1)*H(x,n-2);
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

void list(int N, int D, MatrixXd &order) {
    //Returns the index list used in Slater

    int counter = 0;
    // Two dimensions
    if (D == 2) {
        for(int i=0; i<N; i++) {
            for(int s=i; s<N; s++) {
                int j = s - i;
                order(counter,1) = i;
                order(counter,0) = j;
                counter += 1;
            }
        }
    }

    // Three dimensions
    else if (D == 3) {
        for(int i=0; i<N; i++) {
            for(int j=0; j<N; j++) {
                for(int s=i+j; s<N; s++) {
                    int k = s - i - j;
                    order(counter,0) = i;
                    order(counter,1) = j;
                    order(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }
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


// === CALCULATE DISTANCE MATRIX ===
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


// === ENERGY CALCULATION ===


double energy(const VectorXd &Xa, int D, int O, int k) {
    // Calculate some kinetic energy

    int M = Xa.size();
    int length = M/(2*D);

    MatrixXd A = MatrixXd::Ones(length, length);
    MatrixXd dA = MatrixXd::Zero(length, length);

    if(k<M/2) {
        matrix(Xa.head(M/2), O, D, length, A);
        derivative(Xa.head(M/2), O, D, k, dA);
    }
    else {
        matrix(Xa.tail(M/2), O, D, length, A);
        derivative(Xa.tail(M/2), O, D, k-M/2, dA);
    }


    return (A.inverse()*dA).trace();
}

double energy2(const VectorXd &Xa, int D, int O, VectorXd &diff) {
    // Calculating grad(det(D))

    int M_half = Xa.size()/2;
    int P_half = M_half/D;

    MatrixXd A_up = MatrixXd::Ones(P_half, P_half);
    MatrixXd A_dn = MatrixXd::Ones(P_half, P_half);

    matrix(Xa.head(M_half), O, D, P_half, A_up);
    matrix(Xa.tail(M_half), O, D, P_half, A_dn);

    for(int k = 0; k<M_half; k++) {
        MatrixXd dA = MatrixXd::Zero(P_half, P_half);
        derivative(Xa.head(M_half), O, D, k, dA);

    }
    for(int k = M_half; k<M_half*2; k++) {
        MatrixXd dA = MatrixXd::Zero(P_half, P_half);
        derivative(Xa.tail(M_half), O, D, k, dA);
    }

}
