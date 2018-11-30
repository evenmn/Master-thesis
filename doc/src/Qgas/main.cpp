#include <iostream>
#include <cmath>

#include "vmc.h"
#include "general_tools.h"
#include "test.h"

using namespace std;

int main()
{
    int     P           = 2;               //Number of particles
    int     D           = 2;                //Number of dimensions
    int     N           = P;                //Number of hidden nodes
    int     MC          = int(pow(2,14));        //Number of Monte Carlo cycles
    int     O           = int(orbitals(P,D));    //Number of orbitals
    int     iterations  = 100;             //Number of gradient decent cycles
    int     sampling    = 0;                //Brute force- (0), Hastings- (1) or Gibbs' sampling (2)
    bool    interaction = 1;                //Interaction on if true
    bool    one_body    = 1;                //Calculating onebody density if true
    double  sigma       = 1.0;              //Width of Gaussian distribution
    double  omega       = 1.0;              //Frequency
    double  steplength  = 1.0;              //Steplength for Metropolis
    double  timestep    = 1.0;              //Timestep used in Hastings algorithm
    double  eta         = 0.1; //pow(50*10,-O);     //Learning rate for gradient decent
    double  Diff        = 0.5;              //Diffusion constant


    VMC(P, Diff, D, N, MC, O, iterations, sampling, sigma, omega, steplength, timestep, eta, interaction, one_body);

    //test_without_argument();

    return 0;
}
