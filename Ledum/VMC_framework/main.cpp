#include <iostream>
#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/hydrogenorbital.h"
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
    int numberOfDimensions  = 1;
    int numberOfParticles   = 1;
    int numberOfSteps       = int(1e6);
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.5;          // Variational parameter.
    double beta             = 0.5;          // Variational parameter.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used
    // for equilibration.

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, numberOfParticles, numberOfDimensions));
    system->setWaveFunction             (new SimpleGaussian(system, beta));
    system->setInitialState             (new RandomNormal(system, numberOfDimensions, numberOfParticles));
    system->setOptimizer                (new GradientDescent(system));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps);
    return 0;
}
