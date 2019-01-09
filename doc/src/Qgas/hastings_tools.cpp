#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include "common.h"
#include "eigen3/Eigen/Dense"
#include "wavefunction.h"

#include <random>
#include <cmath>

using namespace Eigen;
using namespace std;

double QForce(const VectorXd &Xa, int i) {

    double QFa = -Xa(i);
    for(int j=0; j<N; j++) {
        QFa += W(i,j)*e(j);
    }

    Slater Slat;
    Jastrow Jast;


    double QF = 0;
    QF += Jast.Jastrow_NQS(v, i, 1);
    QF += Slat.Gauss_ML(Xa, i, 1);
    QF += Slat.SlaterDet(Xa, i, 1);
    QF += Jast.PadeJastrow(i, 1);
    //QF += Slat.Gauss_partly(Xa, k, 1);

    return 2*QF;
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
