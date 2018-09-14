#include <iostream>
#include <random>
#include <cmath>
#include "eigen3/Eigen/Dense"
#include "gradient_descent.h"
#include "test.h"

using namespace std;
using namespace Eigen;


int main()
{
    int     P           = 12;               //Number of particles
    int     D           = 2;               //Number of dimensions
    int     N           = P;              //Number of hidden nodes
    int     MC          = pow(2,16);        //Number of Monte Carlo cycles
    int     iterations  = 3200;             //Number of gradient decent cycles
    int     sampling    = 0;               //Brute force- (0), Hastings- (1) or Gibbs' sampling (2)
    bool    interaction = 1;               //Interaction on if true
    bool    one_body    = 1;               //Calculating onebody density if true
    double  sigma       = 1.0;             //Width of Gaussian distribution
    double  omega       = 1.0;             //Frequency
    double  steplength  = 1.0;             //Steplength for Metropolis
    double  timestep    = 1.0;             //Timestep used in Hastings algorithm
    double  eta         = 0.01;            //Learning rate for gradient decent
    double  Diff        = 0.5;             //Diffusion constant

    GradientDescent(P, Diff, D, N, MC, iterations, sampling, sigma, omega, steplength, timestep, eta, interaction, one_body);

    //test_E_L_calc();

    return 0;
}


