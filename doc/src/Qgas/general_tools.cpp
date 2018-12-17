#include "eigen3/Eigen/Dense"
#include "common.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double binomial(int n, int p) {
    //Binomial coefficients, equal to magic numbers
    return factorial(n+p)/(factorial(n)*factorial(p));
}

double orbitals() {

    int counter = 0;
    while(true) {
        double orb = 2*binomial(counter, D);
        if(int(orb) == P) {
            return int(counter+1);
            break;
        }
        else if(orb > P) {
            std::cout << "Please choose a P such that the orbital is full" << std::endl;
            exit(0);
        }
        counter += 1;
    }
}

string generate_filename(string name, string extension) {

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
