#include <iostream>
#include "wavefunction.h"
#include <random>
#include "eigen3/Eigen/Dense"
#include <cmath>
#include <ctime>
#include <fstream>
#include "hastings_tools.h"
#include "gibbs_tools.h"
#include "test.h"

using namespace Eigen;
using namespace std;

//Mersenne Twister RNG
random_device rd;                   //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                  //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0, 1);
uniform_int_distribution<> hrand(0, 1);

double random_position(){
    return dis(gen);
}

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void GradientDescent(int P, double Diff, int D, int N, int MC, int iterations, int sampling, double sigma, \
                     double omega, double steplength, double timestep, double eta, bool interaction, bool one_body) {

    //Constants
    double psi_ratio = 0;               //ratio of new and old wave function
    double sigma_sqrd = sigma * sigma;
    int M = P*D;
    int M_rand = 0;
    int N_rand = 0;

    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);       //Gaussian distr random number generator
    uniform_int_distribution<> mrand(0, M-1);         //Random number between 0 and M
    uniform_int_distribution<> nrand(0, N-1);         //Random number between 0 and N

    double factor = 0.5;
    double factor_x = 5.0;

    MatrixXd W       = MatrixXd::Random(M, N) * factor;
    VectorXd a       = VectorXd::Random(M)    * factor;
    VectorXd b       = VectorXd::Random(N)    * factor;
    VectorXd X       = VectorXd::Random(M)    * factor_x;
    VectorXd X_new   = VectorXd::Zero(M);
    VectorXd h       = VectorXd::Zero(N);
    VectorXd energies_old = VectorXd::Zero(5);


    VectorXd Xa      = X - a;
    VectorXd v       = b + (W.transpose() * X)/(sigma_sqrd);
    VectorXd X_newa  = VectorXd::Zero(M);
    VectorXd v_new   = VectorXd::Zero(N);

    for(int i=0; i<N; i++) {
        h(i) = hrand(gen);
    }

    WaveFunction Psi;
    Psi.setTrialWF(N, M, D, sampling, sigma_sqrd, omega);

    //Define bins for the one body density measure
    int number_of_bins = 500;
    double max_radius = 5;
    double radial_step = max_radius/number_of_bins;
    double bin_dist[number_of_bins];
    double bins_particles[number_of_bins];

    ofstream ob_file;
    if(one_body) {
        for(int i=0; i<number_of_bins; i++){
            bin_dist[i] = i * radial_step;
            bins_particles[i] = 0;
        }
        ob_file.open ("../../data/ob_density.dat");
    }

    //Open file for writing
    ofstream myfile;
    myfile.open("../../data/energy.txt");

    ofstream myfile1;
    myfile1.open("../../data/local_energies.txt");

    for(int iter=0; iter<iterations; iter++) {
        //averages and energies
        double EL_tot      = 0;          //sum of energies of all states
        double EL_tot_sqrd = 0;          //sum of energies of all states squared
        double E_k         = 0;          //sum of kinetic energies
        double E_ext       = 0;          //sum of potential energy from HO
        double E_int       = 0;          //sum of potential energy from interaction
        double E_k_tot     = 0;
        double E_ext_tot   = 0;
        double E_int_tot   = 0;
        double E = Psi.EL_calc(X, Xa, v, W, D, interaction, E_k, E_ext, E_int);

        VectorXd da_tot           = VectorXd::Zero(M);
        VectorXd daE_tot          = VectorXd::Zero(M);
        VectorXd db_tot           = VectorXd::Zero(N);
        VectorXd dbE_tot          = VectorXd::Zero(N);
        MatrixXd dW_tot           = MatrixXd::Zero(M,N);
        MatrixXd dWE_tot          = MatrixXd::Zero(M,N);

        double accept = 0;
        double tot_dist = 0;

        clock_t start_time = clock();
        for(int i=0; i<MC; i++) {
            X_new = X;              //Setting new matrix equal to old one
            M_rand = mrand(gen);    //Random particle and dimension

            if(sampling==0||sampling==1){
                if(sampling == 0) {
                    //Standard Metropolis
                    X_new(M_rand) = X(M_rand) + (2*random_position() - 1.0)*steplength;
                    X_newa = X_new - a;
                    v_new = b + (W.transpose() * X_new)/sigma_sqrd;
                    psi_ratio = Psi.Psi_value_sqrd(X_newa, v_new)/Psi.Psi_value_sqrd(Xa, v);
                }

                else if(sampling == 1) {
                    //Metropolis-Hastings
                    X_new(M_rand) = X(M_rand) + Diff*QForce(Xa, v, W, sigma_sqrd, M_rand)*timestep + eps_gauss(gen)*sqrt(timestep);
                    X_newa = X_new - a;
                    v_new = b + (W.transpose() * X_new)/sigma_sqrd;
                    psi_ratio = GreenFuncSum(X, X_new, X_newa, Xa, v, W, sigma_sqrd, timestep, D, Diff) * \
                                (Psi.Psi_value_sqrd(X_newa, v_new)/Psi.Psi_value_sqrd(Xa, v));
                }

                if(psi_ratio >= random_position()) {
                    //accept and update
                    accept += 1;
                    X  = X_new;
                    Xa = X_newa;
                    v  = v_new;
                    E  = Psi.EL_calc(X, Xa, v, W, D, interaction, E_k, E_ext, E_int);
                }
            }

            else if(sampling == 2) {
                //Gibbs' sampling
                N_rand = nrand(gen);
                X(M_rand) = x_sampling(a, h, W, sigma_sqrd, M_rand);
                h(N_rand) = h_sampling(v, N_rand);
                Xa = X - a;
                v = b + (W.transpose() * X)/sigma_sqrd;
                E = Psi.EL_calc(X, Xa, v, W, D, interaction, E_k, E_ext, E_int);
            }

            if(one_body && iter == iterations-1) {
                for(int j=0; j<P; j++) {
                    double dist = 0;
                    for(int d=0; d<D; d++) {
                        dist += X(D*j+d)*X(D*j+d);
                    }
                    double r = sqrt(dist);      //Distance from particle to origin
                    for(int k=0; k<number_of_bins; k++) {
                        if(r < bin_dist[k]) {
                            bins_particles[k] += 1;
                            break;
                        }
                    }
                }
            }

            if(iter == iterations - 1) {
                for(int j=0; j<P; j++) {
                    for(int k=0; k<j; k++) {
                        double dist = 0;
                        for(int d=0; d<D; d++)
                            dist += (X(D*j+d) - X(D*k+d))*(X(D*j+d) - X(D*k+d));
                        tot_dist += sqrt(dist);
                    }
                }
            }
            if((iter-iterations-1)%100==0) {
                myfile1 << E << endl;
            }


            VectorXd da = VectorXd::Zero(M);
            VectorXd db = VectorXd::Zero(N);
            MatrixXd dW = MatrixXd::Zero(M,N);
            Psi.Gradient_a(Xa, da);
            Psi.Gradient_b(v, db);
            Psi.Gradient_W(X, v, dW);


            da_tot   += da;
            daE_tot  += E*da;
            db_tot   += db;
            dbE_tot  += E*db;
            dW_tot   += dW;
            dWE_tot  += E*dW;

            EL_tot      += E;
            EL_tot_sqrd += E*E;
            E_k_tot     += E_k;
            E_ext_tot     += E_ext;
            E_int_tot   += E_int;
        }
        clock_t end_time = clock();

        //Calculate <EL> and <EL^2>
        double EL_avg = EL_tot/MC;
        double EL_avg_sqrd = EL_tot_sqrd/MC;
        double variance = (EL_avg_sqrd - EL_avg*EL_avg)/MC;
        double CPU_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;

        cout << "\n--- Iteration " << iter+1 << " ---" << endl;
        cout << "E_L_avg: " << EL_avg << endl;
        cout << "Acceptance ratio: " << accept/MC << endl;
        cout << "Variance " << variance << endl;
        cout << "CPU time: " << CPU_time << "\n" << endl;

        //Printing onebody density to file
        if(iter == iterations - 1) {
            if(one_body){
                //Write to file
                for(int j=0; j<number_of_bins; j++) {
                    ob_file << bins_particles[j] << "\n";
                }
                //Close myfile
                ob_file.close();
            }
            cout << "<E_k>: " << E_k_tot/MC << endl;
            cout << "<E_ext>: " << E_ext_tot/MC << endl;
            cout << "<E_int>: " << E_int_tot/MC << endl;
            cout << "Mean distance: " << tot_dist/(MC*factorial(P-1)) << endl;

            test_energy_convergence(EL_avg, omega, M, interaction);
        }

        //Stop criterion
        for(int j=0; j<4; j++) {
            energies_old(j) = energies_old(j+1);
        }
        energies_old(4) = EL_avg;
        double epsilon = 0.002;

        /*
        if(fabs(EL_avg - energies_old(0))/energies_old(0) < epsilon && \
           fabs(EL_avg - energies_old(1))/energies_old(1) < epsilon && \
           fabs(EL_avg - energies_old(2))/energies_old(2) < epsilon && \
           fabs(EL_avg - energies_old(3))/energies_old(3) < epsilon) {
            cout << "\n Final values" << endl;
            cout << "E_L_avg: " << EL_avg << endl;
            cout << "Acceptance ratio: " << accept/MC << endl;
            cout << "Variance " << variance << endl;
            cout << "CPU time: " << CPU_time << "\n" << endl;
            break;
        }
        */

        //Gradient descent
        a -= 2*eta*(daE_tot - EL_avg*da_tot)/MC;
        b -= 2*eta*(dbE_tot - EL_avg*db_tot)/MC;
        W -= 2*eta*(dWE_tot - EL_avg*dW_tot)/MC;

        //Write to file
        myfile << EL_avg << "\n";
    }
    //Close myfile
    if(myfile.is_open())  myfile.close();
    if(myfile1.is_open()) myfile1.close();

}
