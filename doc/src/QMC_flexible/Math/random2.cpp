#include "random2.h"
#include <random>
#include <iostream>

using namespace std;

//Mersenne Twister RNG
random_device rd;                       //Will be used to obtain a seed for the random number engine
mt19937 gen(rd());                      //Standard mersenne_twister_engine seeded with rd()


double Random2::nextGaussian(double mean, double standardDeviation) {
    normal_distribution<double> rand_gauss(mean, standardDeviation);
    return rand_gauss(gen);
}

int Random2::nextInt(int upperLimit) {
    uniform_int_distribution<> rand_int(0, upperLimit - 1);
    return rand_int(gen);
}

double Random2::nextDouble() {
    uniform_real_distribution<> rand_double(0, 1);
    return rand_double(gen);
}
