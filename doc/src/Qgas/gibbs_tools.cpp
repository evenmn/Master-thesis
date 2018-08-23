#include <iostream>
#include "wavefunction.h"
#include "hastings_tools.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>
#include "gibbs_tools.h"

using namespace std;

random_device RD;                   //Will be used to obtain a seed for the random number engine
mt19937 Gen(RD());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> Dis(0, 1);

double Random_position(){
    return Dis(Gen);
}

double x_sampling(const VectorXd &a, const VectorXd &h, const MatrixXd &W, double sigma_sqrd, int i) {
    double mu = a(i);
    int N = h.size();
    for(int j=0; j<N; j++) {
        mu += h(j)*W(i,j);
    }
    normal_distribution<double> d(mu, sigma_sqrd);
    return d(Gen);

}

double h_sampling(const VectorXd &v, int i){

    double P_h1 = 1/(1 + exp(-v(i)));
    double P_h0 = 1/(1 + exp(v(i)));

    if(Random_position() <= P_h0) {
        return 0;
    } else {
        return 1;
    }
}

