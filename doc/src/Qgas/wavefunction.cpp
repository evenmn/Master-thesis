#include "wavefunction.h"
#include "basis.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include "eigen3/Eigen/Dense"

using namespace std;
using namespace Eigen;

double rij(VectorXd X, int D) {

    double Ep = 0;              // Sum 1/rij
    int P = X.size()/D;

    for(int i=0; i<P; i++) {
        for(int j=0; j<i; j++) {
            double dist = 0;
            for(int d=0; d<D; d++) {
                double diff = X(D*i+d)-X(D*j+d);
                dist += diff*diff;
            }
            Ep += 1/sqrt(dist);
        }
    }
    return Ep;
}

int WaveFunction::setTrialWF(int N, int M, int D, int sampling, double sigma_sqrd, double omega)
{
    m_N          = N;
    m_M          = M;
    m_D          = D;
    m_sampling   = sampling;
    m_sigma_sqrd = sigma_sqrd;
    m_omega_sqrd = omega*omega;
}

double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v)
{
    //Unnormalized wave function

    double Prob = Slater(m_D, Xa, v, m_sigma_sqrd);
    return Prob * Prob;
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

double WaveFunction::EL_calc(const VectorXd X, const VectorXd Xa, const VectorXd v, const MatrixXd W, \
                             int D, int interaction, double &E_kin, double &E_ext, double &E_int) {
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
        diff(i) = energy(Xa, D, i);
    }

    // === ENERGY CALCULATION ===
    // Kinetic energy
    if(m_sampling==2) {
        for(int i=0; i<m_N; i++) {
            E_kin -= 0.5*(double) (Xa.transpose() * W.col(i)) * e_n(i);
            E_kin += 0.5*(double) ((W.col(i)).transpose() * W.col(i)) * e_n(i) * e_p(i);
            for(int j=0; j<m_N; j++) {
                E_kin += 0.25*(double) ((W.col(i)).transpose() * W.col(j)) * e_n(i) * e_n(j);
            }
        }

        E_kin -= 0.5*m_M * m_sigma_sqrd;
        E_kin += 0.25*Xa.transpose() * Xa;
        E_kin = -E_kin/(2 * m_sigma_sqrd * m_sigma_sqrd);
    }

    else {
        /*
        for(int i=0; i<m_N; i++) {
            E_kin += 2*(double) (diff.transpose() * W.col(i)) * e_n(i);
            var += 2*(double) (diff.transpose() * W.col(i)) * e_n(i);
            E_kin -= 2*(double) (Xa.transpose() * W.col(i)) * e_n(i)/m_sigma_sqrd;
            E_kin += (double) ((W.col(i)).transpose() * W.col(i)) * e_n(i)*e_p(i)/m_sigma_sqrd;
            for(int j=0; j<m_N; j++) {
                E_k += (double) ((W.col(i)).transpose() * W.col(j)) * e_n(i) * e_n(j)/m_sigma_sqrd;
            }
        }
        */

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
    if(interaction) E_int = rij(X, D);

    // Harmonic oscillator potential
    E_ext = (double) (X.transpose() * X) * m_omega_sqrd/ 2;

    return E_kin + E_ext + E_int;
}

void WaveFunction::Gradient_a(const VectorXd &Xa, VectorXd &da) {

    if(m_sampling==2) {
        da = 0.5*Xa/m_sigma_sqrd;
    }
    else{
        da = Xa/m_sigma_sqrd;
    }
}

void WaveFunction::Gradient_b(const VectorXd &e, VectorXd &db) {

    if(m_sampling==2) {
        db = 0.5*e;
    }
    else{
        db = e;
    }
}

void WaveFunction::Gradient_W(const VectorXd &X, const VectorXd &e, MatrixXd &dW) {

    if(m_sampling==2) {
        dW = 0.5*X*e.transpose()/m_sigma_sqrd;
    }
    else{
        dW = X*e.transpose()/m_sigma_sqrd;
    }
}
