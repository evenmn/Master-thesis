#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include "common.h"
#include "eigen3/Eigen/Dense"

#include <random>
#include <cmath>

using namespace Eigen;
using namespace std;

double QForce(const VectorXd &Xa, int i) {

    double QF = -Xa(i);
    for(int j=0; j<N; j++) {
        QF += W(i,j)*e(j);
    }
    return QF*(2/(sigma_sqrd));
}

double GreenFuncSum() {
    double GreenSum  = 0;

    for(int i=0; i<P; i++) {
        double GreenFunc = 0;
        for(int j=0; j<D; j++) {
            double QForceOld = QForce(Xa, D*i+j);
            double QForceNew = QForce(X_newa, D*i+j);
            GreenFunc += 0.5*(QForceOld + QForceNew) * (0.5*Diff*dt*(QForceOld - QForceNew)-X_new(D*i+j)+X(D*i+j));
        }
        GreenSum += exp(GreenFunc);
    }

    return GreenSum;
}
