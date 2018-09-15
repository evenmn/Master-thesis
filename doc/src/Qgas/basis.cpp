#include <iostream>
#include <ctime>
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "general_tools.h"

using namespace Eigen;
using namespace std;

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



void matrix(const VectorXd &Xa, int N, int D, MatrixXd &A) {
    //N: Number of fully occupied shells

    int length = binomial(N-1, D);

    MatrixXd order = MatrixXd::Zero(length, D);
    list(N, D, order);


    for(int i=0; i<length; i++) {
        for(int j=0; j<length; j++) {
            for(int k=0; k<D; k++) {
                A(i,j) *= H(Xa(D*i+k), order(j,k));
            }
        }
    }
}


void derivative(const VectorXd &Xa, int N, int D, int k, MatrixXd &dA) {
    //Derivative of A matrix

    int length = binomial(N-1, D);

    MatrixXd order = MatrixXd::Zero(length, D);
    list(N, D, order);

    // Find relevant row
    int row = int(k/D);

    // Find indices of relevant row
    VectorXd a = VectorXd::Zero(D);
    int l = k%D;
    for(int i=0; i<D; i++) {
        a(i) = k-l+i;
    }

    // Find matrix
    for(int i=0; i<length; i++) {
        dA(row, i) = dH(Xa(k), order(i, l));
        for(int j=0; j<D; j++) {
            if(a(j) != k) {
                dA(row, i) *= H(Xa(a(j)), order(i, j));
            }
        }
    }
}


double energy(const VectorXd &Xa, int D, int k) {
    // Calculate some kinetic energy

    int M = Xa.size();
    int N = orbitals(M/D, D);
    int length = M/(2*D);

    MatrixXd A = MatrixXd::Ones(length, length);
    MatrixXd dA = MatrixXd::Zero(length, length);


    if(k<M/2) {
        matrix(Xa.head(M/2), N, D, A);
        derivative(Xa.head(M/2), N, D, k, dA);
    }
    else {
        matrix(Xa.tail(M/2), N, D, A);
        derivative(Xa.tail(M/2), N, D, k-M/2, dA);
    }

    return (A.inverse()*dA).trace();
}


double Slater(int D, const VectorXd &Xa, const VectorXd &v, double sigma_sqrd) {
    // Setting up Slater determinant

    int M = Xa.size();                                // Number of free dimensions
    int P = int(M/D);                                 // Number of particles
    double n_orbitals = orbitals(P,D);                // Number of orbitals

    // Check if the orbitals are fully occupied, otherwise break
    if(fabs(n_orbitals - int(n_orbitals)) > 0.01) {
        cout << "Number of particles needs to be a magic number" << endl;
        exit(0);
    }

    MatrixXd D_up = MatrixXd::Ones(int(P/2),int(P/2));
    MatrixXd D_dn = MatrixXd::Ones(int(P/2),int(P/2));

    matrix(Xa.head(M/2), n_orbitals, D, D_up);
    matrix(Xa.tail(M/2), n_orbitals, D, D_dn);


    return D_up.determinant()*D_dn.determinant()*Gauss_WF(Xa, sigma_sqrd)*Jastrow_NQS(v);
}
