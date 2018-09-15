#include <iostream>


int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double binomial(int n, int p) {
    //Binomial coefficients, equal to magic numbers
    return factorial(n+p)/(factorial(n)*factorial(p));
}

double orbitals(int P, int D) {

    int counter = 0;
    while(true) {
        double orb = 2*binomial(counter, D);
        if(int(orb) == P) {
            return counter+1;
            break;
        }
        else if(orb > P) {
            std::cout << "Please choose a P such that the orbital is full" << std::endl;
            break;
            exit(0);
        }
        counter += 1;
    }
}

/*
char generate_filename(char sampling, int P, int D, int N, int MC, int interaction, double sigma, double omega, double eta, char extention) {
    return char(extention+"_"+sampling+"_"+"P_"+P+)
}
*/
