#include <iostream>
#include <eigen3/Eigen/Dense>
#include "wavefunction.h"
#include "energy.h"
#include "general_tools.h"
#include "basis.h"

using namespace std;
using namespace Eigen;

void test_energy_convergence(double energy, double omega, int M, bool interaction) {

    double error  = 0;
    double tolerance_w_int = 0.1;
    double tolerance_wo_int = 0.05;

    if(interaction) {
        if(M==4) {
            if(fabs(energy - omega * 3) > tolerance_w_int) {
                cout << "Warning: Energy error is larger than the tolerance" << endl;
            }
        }
        else{
            //No analytical values
        }
    }
    else {
        if(fabs(energy - omega * 0.5 * M) > tolerance_wo_int) {
            cout << "Warning: Energy error is larger than the tolerance" << endl;
        }
    }
}

void test_E_L_calc(){
    int M = 2;
    int N = 2;
    int D = 2;
    int norbitals = 1;
    int sampling   = 2;
    double sigma_sqrd = 1;
    double omega = 1.0;
    double omega_sqrd = omega*omega;
    int interaction = 0;

    WaveFunction Psi;

    MatrixXd W       = MatrixXd::Zero(M, N);
    VectorXd a       = VectorXd::Zero(M);
    VectorXd b       = VectorXd::Zero(N);
    VectorXd X       = VectorXd::Zero(M);
    VectorXd v       = VectorXd::Zero(M);
    VectorXd Xa = VectorXd::Zero(M);

    X(0) = 1.0;
    X(1) = 1.0;
    Xa = X - a;
    v = b + (W.transpose() * X)/(sigma_sqrd);


    double E = 0;
    double E_kin = 0;
    double E_ext = 0;
    double E_int = 0;

    //E = Psi.EL_calc(E_kin, E_ext, E_int);
    //cout << "E: " << E << endl;
    //cout << "E_kin: " << E_kin << endl;
    //cout << "E_ext: " << E_ext << endl;
    //cout << "E_int: " << E_int << endl;
}


void test_orbitals() {
    int P = 6;
    int D = 2;

    int orb = orbitals();

    if(orb != 2) {
        cout << "Function 'orbitals' in 'general_tools.cpp' returns wrong answer" << endl;
        exit(0);
    }
}

void test_matrix() {
    int P = 6;          // Particles
    int D = 2;          // Dimensions
    int M = P*D;        // Free dimensions
    int O = 2;          // #orbitals

    MatrixXd A = MatrixXd::Ones(P/2, P/2);
    VectorXd Xa = VectorXd::Random(M);

    Slater Slat;

    Slat.matrix(Xa.head(M/2), H, A);

    MatrixXd B(3,3);
    B << 1, 2*Xa(0), 2*Xa(1),
         1, 2*Xa(2), 2*Xa(3),
         1, 2*Xa(4), 2*Xa(5);

    if(A.isApprox(B) == 0) {
        cout << "Function 'matrix' in 'basis.cpp' returns wrong answer" << endl;
        exit(0);
    }
}

void test_without_argument() {
    test_E_L_calc();
    test_orbitals();
    test_matrix();
}


// === ENERGY CALCULATION ===

/*
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
*/

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
