#include <iostream>
#include <random>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>

#include "eigen3/Eigen/Dense"
#include "wavefunction.h"
#include "hastings_tools.h"
#include "gibbs_tools.h"
#include "general_tools.h"
#include "optimization.h"
#include "test.h"
#include "energy.h"
#include "basis.h"
#include "common.h"

using namespace Eigen;
using namespace std;

//Declare variables connected to common.h
int    M;
int    M_half;
int    P_half;
double sigma_sqrd;

MatrixXd W;
VectorXd X;
VectorXd v;
VectorXd e;
VectorXd e_p;
VectorXd a;
VectorXd b;
VectorXd Xa;
MatrixXd Dist;
MatrixXd Dist_inv;

//Mersenne Twister RNG
random_device rd;                       //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                      //Standard mersenne_twister_engine seeded with rd()
uniform_real_distribution<> dis(0, 1);  //dis gives uniform number between 0 and 1
uniform_int_distribution<> hrand(0, 1); //hrand gives either 0 or 1 randomly

double random_position(){
    return dis(gen);
}

void VMC() {

    //Declare constants
    double psi_ratio = 0;               //ratio of new and old wave function
    sigma_sqrd = sigma * sigma;         //Sigma squared
    M = P*D;                            //Number of free dimensions
    M_half = int(M/2);
    P_half = int(P/2);                  //Number of particles with spin up
    int M_rand = 0;                     //Which M to update
    int N_rand = 0;                     //Which N to update

    //Data path
    string path = "../../data/";        //Path to data folder

    // Make objects
    WaveFunction Psi;                   //Wavefunction object
    Optimization OPT;                   //Optimization object
    Energy ENG;                         //Energy calculation object
    Slater Slat;                        //Slater object

    //Marsenne Twister Random Number Generator
    normal_distribution<double> eps_gauss(0,1);       //Gaussian distr random number generator
    uniform_int_distribution<> mrand(0, M-1);         //Random number between 0 and M
    uniform_int_distribution<> nrand(0, N-1);         //Random number between 0 and N

    double factor = 0.5;                                    //Factor on initial weights
    double factor_x = 5.0;                                  //Factor on initial positions

    W       = MatrixXd::Random(M, N) * factor;     //Initialize W
    a       = VectorXd::Random(M)    * factor;     //Initialize a
    b       = VectorXd::Random(N)    * factor;     //Initialize b
    X       = VectorXd::Random(M)    * factor_x;   //initialize X
    VectorXd X_new   = VectorXd::Zero(M);                   //Declare X_new
    VectorXd h       = VectorXd::Zero(N);                   //Declare h
    e       = VectorXd::Zero(N);                   //Declare e
    VectorXd energies_old = VectorXd::Zero(5);              //Store old energies to check convergence

    Xa      = X - a;                                            //Define useful array Xa = X - a
    v       = b + (W.transpose() * X)/(sigma_sqrd);             //Define useful exponent array
    VectorXd X_newa  = VectorXd::Zero(M);                       //Define updated Xa-array
    VectorXd v_new   = VectorXd::Zero(N);                       //Define updated v-array
    VectorXd e_new   = VectorXd::Zero(N);                       //define updated e-array

    //Set up determinant matrix
    MatrixXd A_up = MatrixXd::Ones(P_half, P_half);         //Slater matrix for spin up
    MatrixXd A_dn = MatrixXd::Ones(P_half, P_half);         //Slater matrix for spin down
    MatrixXd A_up_inv = MatrixXd::Ones(P_half, P_half);     //Slater matrix for spin up inverse
    MatrixXd A_dn_inv = MatrixXd::Ones(P_half, P_half);     //Slater matrix for spin down inverse

    Slat.matrix(Xa.head(M_half), H, A_up);                  //Update Slater matrix for spin up
    Slat.matrix(Xa.tail(M_half), H, A_dn);                  //Update Slater matrix for spin down

    A_up_inv = A_up.inverse();                              //Calculate the inverse
    A_dn_inv = A_dn.inverse();                              //Slater matrices

    // Derivative matrix
    MatrixXd dA_up = MatrixXd::Zero(P,P_half);              //Declare derivatives of the
    MatrixXd dA_dn = MatrixXd::Zero(P,P_half);              //Slater matrices

    ENG.dA_matrix(Xa.head(M_half), dA_up);                  //Update the derivatives of the
    ENG.dA_matrix(Xa.tail(M_half), dA_dn);                  //Slater matrices

    // Distance matrix
    Dist = MatrixXd::Zero(P,P);                             //Declare distance matrix
    Dist_inv = MatrixXd::Zero(P,P);                         //Declare distance matrix
    ENG.rij();                                              //Update distance matrix

    //Update h and e
    for(int i=0; i<N; i++) {
        h(i) = hrand(gen);                                  //Initialize h to be used in Gibbs
        e(i) = 1/(1+exp(-v(i)));                            //Initialize e
    }

    //Define bins for the one body density measure
    int number_of_bins = 500;                               //Number of bins in OB density
    double max_radius = 5;                                  //Maximum bin radius
    double radial_step = max_radius/number_of_bins;         //Bins per radius
    double bin_dist[number_of_bins];                        //Array with start radiuses of the bins
    double bins_particles[number_of_bins];                  //Number of particles in each bin

    ofstream ob_file;
    if(one_body) {
        for(int i=0; i<number_of_bins; i++){
            bin_dist[i] = i * radial_step;                  //Set up bin_dist array
            bins_particles[i] = 0;                          //Initialize bins_particles
        }
        string ob_filename = generate_filename("OB", ".dat");
        ob_file.open (path + ob_filename);                  //Open OB file based on parameters
    }

    //Open file for writing
    string energy_filename = generate_filename("Energy", ".dat");
    string local_filename = generate_filename("Local", ".dat");

    ofstream energy;
    ofstream local;

    energy.open(path + energy_filename);                    //Open energy file based on parameters
    local.open(path + local_filename);                      //Open local energy file based on parameters

    for(int iter=0; iter<iterations; iter++) {
        //averages and energies
        double EL_tot      = 0;          //sum of energies of all states
        double EL_tot_sqrd = 0;          //sum of energies of all states squared
        double E_kin       = 0;          //kinetic energies
        double E_ext       = 0;          //potential energy from HO
        double E_int       = 0;          //potential energy from interaction
        double E_k_tot     = 0;          //sum of kinetic energies
        double E_ext_tot   = 0;          //sum of potential energy from HO
        double E_int_tot   = 0;          //sum of potential energy from interaction
        double accept      = 0;          //Number of accepted moves

        //Initial energy
        double E = ENG.EL_calc(A_up_inv, A_dn_inv, dA_up, dA_dn, E_kin, E_ext, E_int);

        VectorXd da_tot           = VectorXd::Zero(M);      //Declare vector with sum over all da's
        VectorXd daE_tot          = VectorXd::Zero(M);      //Declare vector with sum over all db's
        VectorXd db_tot           = VectorXd::Zero(N);      //Declare vector with sum over all da*E
        VectorXd dbE_tot          = VectorXd::Zero(N);      //Declare vector with sum over all db*E
        MatrixXd dW_tot           = MatrixXd::Zero(M,N);    //Declare matrix with sum over all dW's
        MatrixXd dWE_tot          = MatrixXd::Zero(M,N);    //Declare matrix with sum over all dW*E

        //START MONTE CARLO INTEGRATION
        clock_t start_time = clock();
        for(int i=0; i<MC; i++) {

            X_new = X;              //Setting new matrix equal to old one
            M_rand = mrand(gen);    //Random particle and dimension

            if(sampling==0||sampling==1){
                //Metropolis sampling
                if(sampling == 0) {
                    //Standard Metropolis
                    X_new(M_rand) = X(M_rand) + (2*random_position() - 1.0)*steplength;         //Update position
                    X_newa = X_new - a;                                                         //Update X - a
                    v_new = b + (W.transpose() * X_new)/sigma_sqrd;                             //Update v
                    //for (int i=0; i<N; i++) e_new(i) = 1/(1+exp(-v_new(i)));                    //Update e
                    psi_ratio = Psi.Psi_value_sqrd(X_newa, v_new)/Psi.Psi_value_sqrd(Xa, v);    //Calculate ratio
                }

                else if(sampling == 1) {
                    //Metropolis-Hastings
                    X_new(M_rand) = X(M_rand) + Diff*QForce(Xa, v, W, sigma_sqrd, M_rand) * \
                                    timestep + eps_gauss(gen)*sqrt(timestep);                   //Update position
                    X_newa = X_new - a;                                                         //Update X - a
                    v_new = b + (W.transpose() * X_new)/sigma_sqrd;                             //Update v
                    //for (int i=0; i<N; i++) e_new(i) = 1/(1+exp(-v_new(i)));                    //Update e
                    psi_ratio = GreenFuncSum(X, X_new, X_newa, Xa, v, W, sigma_sqrd, timestep, D, Diff) * \
                                (Psi.Psi_value_sqrd(X_newa, v_new)/Psi.Psi_value_sqrd(Xa, v));  //Calculate ratio
                }

                if(psi_ratio >= random_position()) {
                    //accept and update
                    accept += 1;                        //Add 1 to number of accepted moves
                    X  = X_new;                         //Set new position to actual position
                    Xa = X_newa;                        //Set new Xa to actual Xa
                    v  = v_new;                         //Set new v to actual v
                    //for (int i=0; i<N; i++) e(i) = 1/(1+exp(-v(i)));                    //Update e
                    //e  = e_new;                         //Set new e to actual e

                    int row = int(M_rand/D);            //Which row in the Slater matrix to update
                    ENG.rij_cross(row);                 //Update distance matrix

                    // Find indices of relevant row
                    VectorXd c = VectorXd::Zero(D);     //Vector with all dimensions
                    int l = M_rand%D;                   //With particle to update position to
                    for(int i=0; i<D; i++) {
                        c(i) = M_rand-l+i;              //Update vector with all dimensions in correct order
                    }

                    double R = 0;                       //To be used when calculate inverse iteratively
                    if(row < P_half){
                        Slat.A_rows(Xa.head(M_half), H, row, A_up);     //Update Slater matrix

                        A_up_inv = A_up.inverse();
                        /*
                        for(int j=0; j<P_half; j++) {
                            R += A_up(row, j)*A_up_inv(j, row);
                        }
                        for(int j=0; j<P_half; j++) {
                            A_up_inv(j, row) /= R;              //Update inverse of Slater matrix by iterations
                        }
                        */
                        for(int i=0; i<D; i++) {
                            ENG.dA_row(Xa.head(M_half), int(c(i)), dA_up);         //Update derivative of Slater matrix
                        }

                    }
                    else {
                        row -= P_half;
                        Slat.A_rows(Xa.tail(M_half), H, row, A_dn); //Update Slater matrix

                        A_dn_inv = A_dn.inverse();

                        /*
                        for(int j=0; j<P_half; j++) {
                            R += A_dn(row, j)*A_dn_inv(j,row);
                        }
                        for(int j=0; j<P_half; j++) {
                            A_dn_inv(j,row) /=R;                        //Update inverse of Slater matrix by iterations
                        }
                        */
                        for(int i=0; i<D; i++) {
                            ENG.dA_row(Xa.tail(M_half), int(c(i)-M_half), dA_dn);     //Update derivative of Slater matrix
                        }
                    }

                    E  = ENG.EL_calc(A_up_inv, A_dn_inv, dA_up, dA_dn, E_kin, E_ext, E_int);    //Update energy
                }
            }

            else if(sampling == 2) {
                //Gibbs' sampling
                N_rand = nrand(gen);                                    //Random hidden node
                X(M_rand) = x_sampling(a, h, W, sigma_sqrd, M_rand);    //Update position of random free dimension
                h(N_rand) = h_sampling(v, N_rand);                      //Update random hidden node
                Xa = X - a;                                             //Update Xa
                v = b + (W.transpose() * X)/sigma_sqrd;                 //Update v

                //Additional stuff
                int row = int(M_rand/D);                                //Which row in the Slater matrix to update

                // Find indices of relevant row
                VectorXd c = VectorXd::Zero(D);     //Vector with all dimensions
                int l = M_rand%D;                   //With particle to update position to
                for(int i=0; i<D; i++) {
                    c(i) = M_rand-l+i;              //Update vector with all dimensions in correct order
                }

                double R = 0;                       //To be used when calculate inverse iteratively
                if(row < P/2){
                    Slat.A_rows(Xa.head(M_half), H, row, A_up);     //Update Slater matrix

                    A_up_inv = A_up.inverse();
                    /*
                    for(int j=0; j<P_half; j++) {
                        R += A_up(row, j)*A_up_inv(j, row);
                    }
                    for(int j=0; j<P_half; j++) {
                        A_up_inv(j, row) /= R;              //Update inverse of Slater matrix by iterations
                    }
                    */
                    for(int i=0; i<D; i++) {
                        ENG.dA_row(Xa.head(M_half), int(c(i)), dA_up);         //Update derivative of Slater matrix
                    }

                }
                else {
                    row -= P_half;
                    Slat.A_rows(Xa.tail(M_half), H, row, A_dn); //Update Slater matrix

                    A_dn_inv = A_dn.inverse();
                    /*
                    for(int j=0; j<P_half; j++) {
                        R += A_dn(row, j)*A_dn_inv(j,row);
                    }
                    for(int j=0; j<P_half; j++) {
                        A_dn_inv(j,row) /=R;                        //Update inverse of Slater matrix by iterations
                    }
                    */
                    for(int i=0; i<D; i++) {
                        ENG.dA_row(Xa.tail(M_half), int(c(i)-M_half), dA_dn);     //Update derivative of Slater matrix
                    }
                }

                E  = ENG.EL_calc(A_up_inv, A_dn_inv, dA_up, dA_dn, E_kin, E_ext, E_int);    //Update energy
            }

            //CALCULATE ONEBODY DENSITIES
            if(one_body && iter == iterations-1) {
                for(int j=0; j<P; j++) {
                    double dist = 0;
                    for(int d=0; d<D; d++) {
                        dist += X(D*j+d)*X(D*j+d);
                    }
                    double r = sqrt(dist);      //Distance from particle j to origin
                    for(int k=0; k<number_of_bins; k++) {
                        if(r < bin_dist[k]) {
                            bins_particles[k] += 1;
                            break;
                        }
                    }
                }
            }

            //WRITE LOCAL ENERGY TO FILE FOR BLOCKING
            if((iter-iterations-1)%100==0) {
                local << E << endl;
            }

            //UPDATE GRADIENTS
            VectorXd da = VectorXd::Zero(M);
            VectorXd db = VectorXd::Zero(N);
            MatrixXd dW = MatrixXd::Zero(M,N);
            OPT.Gradient_a(Xa, da);
            OPT.Gradient_b(e, db);
            OPT.Gradient_W(X, e, dW);


            da_tot   += da;
            daE_tot  += E*da;
            db_tot   += db;
            dbE_tot  += E*db;
            dW_tot   += dW;
            dWE_tot  += E*dW;

            EL_tot      += E;
            EL_tot_sqrd += E*E;
            E_k_tot     += E_kin;
            E_ext_tot   += E_ext;
            E_int_tot   += E_int;
        }   // End Monte Carlo
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
            cout << "Mean distance: " << Dist.sum()/(MC*factorial(P-1)) << endl;

            test_energy_convergence(EL_avg, omega, M, interaction);
        }

        //Stop criterion
        /*
        for(int j=0; j<4; j++) {
            energies_old(j) = energies_old(j+1);
        }
        energies_old(4) = EL_avg;
        double epsilon = 0.002;

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

        // OPTIMIZATION
        // Declare optimization arrays
        VectorXd opt_a = VectorXd::Zero(M);
        VectorXd opt_b = VectorXd::Zero(N);
        MatrixXd opt_W = MatrixXd::Zero(M,N);

        if(optimization==0) {
            // Gradient Descent
            OPT.GD_a(EL_avg, daE_tot, da_tot, opt_a);
            OPT.GD_b(EL_avg, dbE_tot, db_tot, opt_b);
            OPT.GD_W(EL_avg, dWE_tot, dW_tot, opt_W);
        }
        else if(optimization==1) {
            // ADAM
            VectorXd ma = VectorXd::Zero(M);
            VectorXd va = VectorXd::Zero(M);
            VectorXd mb = VectorXd::Zero(N);
            VectorXd vb = VectorXd::Zero(N);
            MatrixXd mW = MatrixXd::Zero(M,N);
            MatrixXd vW = MatrixXd::Zero(M,N);
            double b1 = 0.9;
            double b2 = 0.99;
            OPT.ADAM_a(iter, ma, va, b1, b2, EL_avg, daE_tot, da_tot, opt_a);
            OPT.ADAM_b(iter, mb, vb, b1, b2, EL_avg, dbE_tot, db_tot, opt_b);
            OPT.ADAM_W(iter, mW, vW, b1, b2, EL_avg, dWE_tot, dW_tot, opt_W);
        }

        // Optimization
        a -= opt_a;
        b -= opt_b;
        W -= opt_W;

        //Write to file
        energy << EL_avg << "\n";
    }   // End iterations
    //Close myfile
    if(energy.is_open())  energy.close();
    if(local.is_open()) local.close();

}
