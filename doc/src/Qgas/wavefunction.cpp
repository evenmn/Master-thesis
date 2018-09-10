#include "wavefunction.h"
#include "basis.h"
#include <cmath>
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
                dist += (X(D*i+d)-X(D*j+d))*(X(D*i+d)-X(D*j+d));
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

    double Prob = Slater(m_M/m_D, m_D, Xa, v, m_sigma_sqrd);
    return Prob * Prob;
}

void Deter(const VectorXd &Xa, VectorXd &diff) {
    // Determinant dependent part
    MatrixXd D_up = MatrixXd::Zero(int(3),int(3));
    MatrixXd D_dn = MatrixXd::Zero(int(3),int(3));

    D_up << H(0,0), H(Xa[0],1), H(Xa[1],1),
            H(0,0), H(Xa[2],1), H(Xa[3],1),
            H(0,0), H(Xa[4],1), H(Xa[5],1);

    D_dn << H(0,0), H(Xa[6],1), H(Xa[7],1),
            H(0,0), H(Xa[8],1), H(Xa[9],1),
            H(0,0), H(Xa[10],1), H(Xa[11],1);


    VectorXd X_up = VectorXd::Zero(12);
    VectorXd X_dn = VectorXd::Zero(12);
    for(int i=0; i<6; i++) {
        X_up(i) = X_up(i+6) = Xa(i);
        X_dn(i) = X_dn(i+6) = Xa(i+6);
    }

    for(int i=0; i<6; i++) {
        if(i % 2==0){
            diff(i) = 4*(X_up(i+3) - X_up(i+5))/D_up.determinant();
            diff(i+6) = 4*(X_dn(i+3) - X_dn(i+5))/D_dn.determinant();
        }
        else {
            diff(i) = 4*(X_up(i+3) - X_up(i+1))/D_up.determinant();
            diff(i+6) = 4*(X_dn(i+3) - X_dn(i+1))/D_dn.determinant();
        }
    }
}

double WaveFunction::EL_calc(const VectorXd X, const VectorXd Xa, const VectorXd v, const MatrixXd W, \
                             int D, int interaction, double &E_k, double &E_ext, double &E_int) {
    // Local energy calculations

    double E = 0;
    E_k = 0;
    E_ext = 0;
    E_int = 0;
    double E_knew = 0;
    double E_pnew = 0;
    double E_intnew = 0;
    double number = 0;

    VectorXd e = VectorXd::Zero(m_N);
    VectorXd eNominator = VectorXd::Zero(m_N);
    for(int i=0; i<m_N; i++) {
        double expi = exp(-v(i));
        eNominator(i) = expi;
        e(i) = 1/(1 + expi);
    }

    VectorXd diff = VectorXd::Zero(m_M);
    Deter(Xa, diff);

    // Kinetic energy
    if(m_sampling==2) {
        for(int i=0; i<m_N; i++) {
            E_knew -= 0.5*(double) (Xa.transpose() * W.col(i)) * e(i);
            E_knew += 0.5*(double) ((W.col(i)).transpose() * W.col(i)) * eNominator(i) * e(i) * e(i);
            for(int j=0; j<m_N; j++) {
                E_knew += 0.25*(double) ((W.col(i)).transpose() * W.col(j)) * e(i) * e(j);
            }
        }

        E_knew -= 0.5*m_M * m_sigma_sqrd;
        E_knew += 0.25*Xa.transpose() * Xa;
        E_knew = -E_knew/(2 * m_sigma_sqrd * m_sigma_sqrd);
        E_k += E_knew;
    }

    else {
        for(int i=0; i<m_N; i++) {
            E_knew += 2*(double) (diff.transpose() * W.col(i)) * e(i);
            E_knew -= 2*(double) (Xa.transpose() * W.col(i)) * e(i)/m_sigma_sqrd;
            E_knew += (double) ((W.col(i)).transpose() * W.col(i)) * eNominator(i) * e(i)*e(i)/m_sigma_sqrd;
            for(int j=0; j<m_N; j++) {
                E_knew += (double) ((W.col(i)).transpose() * W.col(j)) * e(i) * e(j)/m_sigma_sqrd;
            }
        }

        E_knew -= m_M;
        E_knew += (double) (Xa.transpose() * Xa)/m_sigma_sqrd;
        E_knew -= 2*(double) (diff.transpose()*Xa);
        //E_knew += (double) (diff.transpose()*diff);
        E_knew = -E_knew/(2 * m_sigma_sqrd);
        //E_k += E_knew;
    }


    // Interaction energy
    if(interaction) E_intnew += rij(X, D);
    E_int += E_intnew;

    // Harmonic oscillator potential
    E_pnew += (double) (X.transpose() * X) * m_omega_sqrd/ 2;
    E_ext += E_pnew;

    E = E_knew + E_pnew + E_intnew;
    E_k = E_knew;
    E_ext = E_pnew;
    E_int = E_intnew;
    return E;
}

void WaveFunction::Gradient_a(const VectorXd &Xa, VectorXd &da) {

    VectorXd diff = VectorXd::Zero(m_M);
    Deter(Xa, diff);

    if(m_sampling==2) {
        da = 0.5*Xa/m_sigma_sqrd-diff;
    }
    else{
        da = Xa/m_sigma_sqrd-diff;
    }

}

void WaveFunction::Gradient_b(const VectorXd &v, VectorXd &db) {

    if(m_sampling==2) {
        for(int i=0; i<m_N; i++)
            db(i) = 0.5/(1 + exp(-v(i)));
    }
    else{
        for(int i=0; i<m_N; i++)
            db(i) = 1/(1 + exp(-v(i)));
    }
}

void WaveFunction::Gradient_W(const VectorXd &X, const VectorXd &v, MatrixXd &dW) {

    if(m_sampling==2) {
        for(int i=0; i<m_N; i++) {
            for(int j=0; j<m_M; j++) {
                dW(j,i) = 0.5*X(j)/(m_sigma_sqrd*(1 + exp(-v(i))));
            }
        }
    }
    else{
        for(int i=0; i<m_N; i++) {
            for(int j=0; j<m_M; j++) {
                dW(j,i) = X(j)/(m_sigma_sqrd*(1 + exp(-v(i))));
            }
        }
    }
}
