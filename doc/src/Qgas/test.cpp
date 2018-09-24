#include <iostream>
#include <eigen3/Eigen/Dense>
#include "wavefunction.h"
#include "basis.h"
#include "general_tools.h"

using namespace std;
using namespace Eigen;

void test_without_argument() {
    test_E_L_calc();
    test_orbitals();
    test_matrix();
}

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
    Psi.setTrialWF(N, M, D, norbitals, sampling, sigma_sqrd, omega);

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

    //E = Psi.EL_calc(X, Xa, v, W, interaction, E_kin, E_ext, E_int);
    //cout << "E: " << E << endl;
    //cout << "E_kin: " << E_kin << endl;
    //cout << "E_ext: " << E_ext << endl;
    //cout << "E_int: " << E_int << endl;
}

/*Should implement tests which test diff for 2 particles and that the matrices are set up correctly for 6 particles
 */

void test_orbitals() {
    int P = 6;
    int D = 2;

    int orb = orbitals(P,D);

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


    matrix(Xa.head(M/2), O, D, P/2, A);

    MatrixXd B(3,3);
    B << 1, 2*Xa(0), 2*Xa(1),
         1, 2*Xa(2), 2*Xa(3),
         1, 2*Xa(4), 2*Xa(5);

    if(A.isApprox(B) == 0) {
        cout << "Function 'matrix' in 'basis.cpp' returns wrong answer" << endl;
        exit(0);
    }
}
