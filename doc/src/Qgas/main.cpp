#include <iostream>
#include <cmath>

#include "vmc.h"
#include "general_tools.h"
#include "test.h"
#include "common.h"

using namespace std;

int     P;
int     D;
int     N;
int     MC;
int     O;
int     iterations;
int     sampling;
int     optimization;
bool    interaction;
bool    one_body;
double  sigma;
double  omega;
double  dx;
double  dt;
double  eta;
double  Diff;

int main()
{
    P           = 2;                //Number of particles
    D           = 2;                //Number of dimensions
    N           = P;                //Number of hidden nodes
    MC          = int(pow(2,18));   //Number of Monte Carlo cycles
    O           = orbitals();       //Number of orbitals
    iterations  = 1000;             //Number of gradient decent iterations
    sampling    = 1;                //Brute force- (0), Hastings- (1) or Gibbs' sampling (2)
    optimization= 0;                //Gradient Descent (0), ADAM (1)
    interaction = 1;                //Interaction on if true
    one_body    = 1;                //Calculating onebody density if true
    sigma       = 1.0;              //Width of Gaussian distribution
    omega       = 1.0;              //Frequency
    dx          = 1.0;              //Steplength for Metropolis
    dt          = 0.1;              //Timestep used in Hastings algorithm
    eta         = 5*pow(10,-O);     //Learning rate for optimization
    Diff        = 0.5;              //Diffusion constant

    VMC();

    //test_without_argument();

    return 0;
}
