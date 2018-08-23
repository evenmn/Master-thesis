#include <iostream>
#include <ctime>
#include <cmath>
#include "eigen3/Eigen/Dense"

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


int magic_numbers(int n) {
    // Magic numbers of a D-dimensional quantum dot.
    // Returns total number of particles, N, in n fully
    // occupied orbitals. Starts at 0, following the
    // arithmetic series
    //
    //    N(n) = S * binomial(n+D, D)
    //
    //with S as the number of spin configurations.

    if(n < 0) {
        return 0;
    }
    else {
        return (n+1)*(n+2);
    }
}

double magic_numbers_inverse(int sum) {
    // Given a magic number, 'sum', this function
    // returns number of orbitals. To be used when
    // checking if orbitals are fully occupied.

    return (-3 + sqrt(1 + 4*sum))/2;
}

double Slater(int P, int D, VectorXd Xa, VectorXd v, double sigma_sqrd) {
    // Setting up Slater determinant

    double n_orbitals = magic_numbers_inverse(P);   // Number of orbitals given N

    // Check if the orbitals are fully occupied, otherwise break
    if(fabs(n_orbitals - int(n_orbitals)) > 0.01) {
        cout << "Number of particles needs to be a magic number" << endl;
        exit(0);
    }
    /*
    MatrixXd D_up = MatrixXd::Zero(int(P/2),int(P/2));
    MatrixXd D_dn = MatrixXd::Zero(int(P/2),int(P/2));

    D_up << H(0,0), H(Xa[0],1), H(Xa[1],1),
            H(0,0), H(Xa[4],1), H(Xa[5],1),
            H(0,0), H(Xa[8],1), H(Xa[9],1);

    D_dn << H(0,0), H(Xa[2],1), H(Xa[3],1),
            H(0,0), H(Xa[6],1), H(Xa[7],1),
            H(0,0), H(Xa[10],1), H(Xa[11],1);


    double n_something = 2*(n_orbitals + 1);

    VectorXd levels = VectorXd::Zero(n_something);
    for(int i = 0; i<(n_orbitals+1); i++) {
        levels(2*i) = i;
        levels(2*i+1) = i;
    }

    MatrixXd indices = MatrixXd::Zero(P, D);
    int counter = 0;
    for(int i = 0; i<n_something; i++) {
        for(int j = 0; j<levels(i)+1; j++) {
            indices(counter,0) = j;
            indices(counter,1) = levels(i) - j;
            counter += 1;
        }
    }

    MatrixXd A = MatrixXd::Zero(P,P);
    for(int i = 0; i<P; i++) {
        for(int j = 0; j<P; j++) {
            A(i,j) = hermite(Xa(2*i), indices(j,0)) * hermite(Xa(2*i+1), indices(j,1));
        }
    }
    cout << A << endl;
    return A.determinant() * NQS_WF(Xa, v, sigma);
    */
    return Gauss_WF(Xa, sigma_sqrd)*Jastrow_NQS(v);//D_up.determinant()* D_dn.determinant()*Gauss_WF(Xa, sigma_sqrd)*Jastrow_NQS(v);
}

/*
int main() {
    int P = 6;      // # Particles
    int D = 2;      // # Dimensions
    int M = P*D;    // # Visible nodes
    int N = 2;      // # Hidden nodes
    int sigma = 1;  // # Variance

    VectorXd Xa = VectorXd::Zero(M);
    VectorXd v = VectorXd::Zero(N);

    Xa << 0.3, 0.5, 0.1, 0.4, 0.7, 0.2, 0.6, 0.9, 0.2, 0.3, 0.8, 0.8;
    v << 1, 1;

    clock_t start_time = clock();
    cout << Slater(P, D, Xa, v, sigma) << endl;
    clock_t end_time = clock();

    double CPU_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    cout << "CPU time: " << CPU_time << "\n" << endl;

    return 0;
}
*/
