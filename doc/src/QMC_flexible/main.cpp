#include <iostream>
#include <Eigen/Dense>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/hydrogenorbital.h"
#include "WaveFunctions/gausspadejastrow.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/atomicnucleus.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/randomnormal.h"
#include "Optimization/optimization.h"
#include "Optimization/gradientdescent.h"
#include "Math/random.h"

using namespace std;


int main() {
    int numberOfDimensions  = 3;
    int numberOfParticles   = 2;
    int numberOfSteps       = int(1e6);
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 1.0;          // Variational parameter.
    double beta             = 1.0;          // Variational parameter.
    double stepLength       = 0.1;          // Metropolis step length.
    bool interaction        = false;
    double equilibration    = 0.1;          // Amount of the total steps used
    // for equilibration.
    
    Eigen::MatrixXd Gamma   = Eigen::MatrixXd::Ones(numberOfParticles, numberOfParticles);

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, numberOfParticles, numberOfDimensions, interaction));
    system->setWaveFunction             (new SimpleGaussian(system, alpha));
    system->setInitialState             (new RandomNormal(system, numberOfDimensions, numberOfParticles));
    system->setOptimizer                (new GradientDescent(system));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps);
    return 0;
}
