#include "wavefunction.h"
#include "general_tools.h"
#include "basis.h"
#include "eigen3/Eigen/Dense"
#include "common.h"

#include <cstdarg>   // Variadic templates
#include <cmath>     // Math
#include <ctime>     // Time
#include <iostream>  // Function

using namespace std;
using namespace Eigen;

// === SLATER ===

double Slater::A_elements(const VectorXd &Xa, double f(double, int), int i, int j) {
    // Updates an element in A-matrix

    MatrixXd order = MatrixXd::Zero(P_half, D);
    list(O, D, order);

    double element = 1;
    for(int k=0; k<D; k++) {
        element *= f(sqrt(omega) * Xa(D*i+k), int(order(j,k)));
    }

    return element;
}

void Slater::A_rows(const VectorXd &Xa, double f(double, int), int i, MatrixXd &A) {
    // Updates a row in A-matrix

    for(int j=0; j<P_half; j++) {
        A(i,j) = A_elements(Xa, f, i, j);
    }
}

void Slater::matrix(const VectorXd &Xa, double f(double, int), MatrixXd &A) {
    // Update the entire matrix

    for(int j=0; j<P_half; j++) {
        A_rows(Xa, f, j, A);
    }
}

double Slater::Gauss(const VectorXd &X, double alpha, int k, int type) {
    // Gaussian with a variational parameter

    if(type == 0) {
        return exp(-double(omega * X.transpose() * X) * alpha);
    }
    else if(type == 1) {
        // First derivative
        return -double(omega * X(k))/alpha;
    }
    else if(type == 2) {
        // Second derivative
        return -omega/alpha;
    }
}


double Slater::Gauss_ML(const VectorXd &Xa, int k, int type) {
    // Biased Gaussian
    if(type == 0) {
        return exp(-double(omega * Xa.transpose() * Xa)/(2 * sigma_sqrd));
    }
    else if(type == 1) {
        // First derivative
        return -double(omega * Xa(k))/sigma_sqrd;
    }
    else if(type == 2) {
        // Second derivative
        return -M/sigma_sqrd;
    }
}


double Slater::Gauss_partly(const VectorXd &X, int k, int type) {
    // Part to go from restricted to partly restricted
    if(type == 0) {
        return exp(-double(X.transpose()*C*X)/(2*sigma_sqrd));
    }
    else if(type == 1) {
        //First derivative
        double result = 0;
        result += C(k,k) * X(k);
        for(int i=0; i<M; i++) {
            result += C(k,i) * X(i);
        }
        return -result/(2*sigma_sqrd);
    }
    else if(type == 2) {
        //Second derivative
        double result = 0;
        for(int i=0; i<M; i++) {
            result += C(i,i);
        }
        return -result/sigma_sqrd;
        //return -(C.diagonal()).sum()/sigma_sqrd;
    }
}


double Slater::SlaterDet(const VectorXd &Xa, int k, int type) {
    // Setting up Slater determinant

    if(type == 0) {
        MatrixXd D_up = MatrixXd::Ones(P_half,P_half);
        MatrixXd D_dn = MatrixXd::Ones(P_half,P_half);

        matrix(Xa.head(M_half), H, D_up);
        matrix(Xa.tail(M_half), H, D_dn);

        return D_up.determinant()*D_dn.determinant();
    }
    else if(type == 1) {
        // First derivative
        return diff(k);
    }
    else if(type == 2) {
        // Second derivative
        return -double(diff.transpose() * diff);
    }
}


// === JASTROW ===

double Jastrow::Jastrow_NQS(const VectorXd &v, int k, int type) {
    //Neural Quantum State Wavefunction (NQS-WF)

    if(type == 0) {
        double prod = 1;
        for(int i=0; i < N; i++) {
            prod *= (1 + exp(v(i)));
        }
        return prod;
    }
    else if(type == 1) {
        // First derivative
        return double(W.row(k) * e)/sigma_sqrd;
    }
    else if(type == 2) {
        // Second derivative
        return (W.cwiseAbs2()*e_p.cwiseProduct(e)).sum();
    }
}


double Jastrow::PadeJastrow(int k, int type) {
    //Pade-Jastrow factor
    MatrixXd A = MatrixXd::Ones(P, P);
    double B = 1;
    if(type == 0) {
        double sum = 0;
        for(int i=0; i<P; i++) {
            for(int j=0; j<i; j++) {
                double rij = Dist(i,j);
                double f = 1/(1 + B*rij);
                sum += A(i,j)*f*rij;
            }
        }
        return exp(sum);
    }
    else if(type == 1) {
        // First derivative
        double result = 0;
        int k_p = int(k/D);     //Particle
        int k_d = k%D;          //Dimension
        for(int j=0; j<k_p; j++) {
            double ximxj = X(k)-X(D*j+k_d);
            double rij = Dist(k_p,j);
            double f = 1/(1 + B*rij);
            result += A(k_p,j)*f*f*ximxj/rij;
        }
        return result;
    }
    else if(type == 2) {
        // Second derivative
        double result = 0;
        for(int i=0; i<M; i++) {
            int i_p = int(i/D);     //Particle
            int i_d = i%D;          //Dimension
            for(int j=0; j<i_p; j++) {
                double ximxj = (X(i)-X(D*j+i_d))*(X(i)-X(D*j+i_d));
                double rij = Dist(i_p,j);
                double f = 1/(1 + B*rij);
                result += (A(i_p,j)*f*f/rij)*(1 - ximxj/(rij*rij) - 2*B*ximxj*f);
            }
        }
        return result;
    }
}


double Jastrow::CartPadeJastrow(int k, int type) {
    //Cartesian Pade-Jastrow
    MatrixXd A = MatrixXd::Ones(P, P);
    double B = 1;
    if(type == 0) {
        double sum = 0;
        for(int i=0; i<P; i++) {
            for(int j=0; j<i; j++) {
                for(int d=0; d<D; d++) {
                    double dist = (X(D*i+d)-X(D*j+d));
                    double f = 1/(1+B*dist);
                    sum += A(i,j)*dist*f;
                }
            }
        }
        return exp(sum);
    }
    else if(type == 1) {
        // First derivative
        double sum = 0;
        int k_p = int(k/D);     //Particle
        int k_d = k%D;          //Dimension
        for(int j=0; j<k_p; j++) {
            double dist = (X(k)-X(D*j+k_d));
            double f = 1/(1+B*dist);
            sum += A(k_p,j)*f*f;
        }
        return sum;
    }
    else if(type == 2) {
        // Second derivative
        double sum = 0;
        for(int k=0; k<M; k++) {
            int k_p = int(k/D);     //Particle
            int k_d = k%D;          //Dimension
            for(int j=0; j<k_p; j++) {
                double dist = (X(k)-X(D*j+k_d));
                double f = 1/(1+B*dist);
                sum -= 2*A(k_p,j)*B*f*f*f;
            }
        }
        return sum;
    }
}


double WaveFunction::Psi_value(const VectorXd &Xa, const VectorXd &v) {
    // Setting up total wavefunction

    Slater Slat;
    Jastrow Jast;

    double slater = 0;
    double jastrow = 0;

    int system = 0;
    if(system == 0) {
        // Hermite functions and NQS

        slater = Slat.SlaterDet(Xa, 0, 0) * Slat.Gauss_ML(Xa, 0, 0); // * Slat.Gauss_partly(Xa, 0, 0);
        jastrow = Jast.PadeJastrow(0, 0) * Jast.Jastrow_NQS(v, 0, 0);
    }

    return slater * jastrow;
}

double WaveFunction::Psi_value_sqrd(const VectorXd &Xa, const VectorXd &v) {
    //Unnormalized wave function squared

    double Prob = Psi_value(Xa, v);
    return Prob * Prob;
}

double WaveFunction::Psi_ratio() {
    //Ratio between new and old probability distribution

    return Psi_value_sqrd(X_newa, v_new)/Psi_value_sqrd(Xa, v);
}
