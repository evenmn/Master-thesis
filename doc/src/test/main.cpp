#include <iostream>
#include <ctime>
#include "eigen3/Eigen/Dense"
#include <string>

using namespace std;
using namespace Eigen;


string generate_filename(int sampling, int P, int D, int N, int MC, int interaction, double sigma, double omega, double eta, string name, string extension) {

    // Sampling method
    string str1;
    if(sampling == 0) {
        str1 = "_BF";
    }
    else if(sampling == 1) {
        str1 = "_H";
    }
    else if(sampling == 2) {
        str1 = "_G";
    }

    // Interaction
    string str2;
    if(interaction == 0) {
        str2 = "_NoInt";
    }
    else if(interaction == 1) {
        str2 = "_Int";
    }

    string str3 = "_P_" + to_string(P);
    string str4 = "_D_" + to_string(D);
    string str5 = "_N_" + to_string(N);
    string str6 = "_MC_" + to_string(MC);

    // Sigma
    int sigma_int = int(sigma);
    int sigma_dec = int(fabs(sigma - sigma_int)*100);
    string str7 = "_sigma_" + to_string(sigma_int) + "p" + to_string(sigma_dec);

    // Omega
    int omega_int = int(omega);
    int omega_dec = int(fabs(omega - omega_int)*100);
    string str8 = "_omega_" + to_string(omega_int) + "p" + to_string(omega_dec);

    // Eta
    int eta_int = int(eta);
    int eta_dec = int(fabs(eta - eta_int)*1000);
    string str9 = "_eta_" + to_string(eta_int) + "p" + to_string(eta_dec);

    return name + str1 + str2 + str3 + str4 + str5 + str6 + str7 + str8 + str9 + extension;
}

int main()
{
    int sampling = 1;
    int P = 2;
    int D = 2;
    int N = 2;
    int MC = 1000;
    int interaction = 1;
    double sigma = 1.0;
    double omega = 1.0;
    double eta = 0.1;
    string name = "energy";
    string extension = ".dat";

    cout << generate_filename(sampling, P, D, N, MC, interaction, sigma, omega, eta, name, extension) << endl;

    return 0;
}
