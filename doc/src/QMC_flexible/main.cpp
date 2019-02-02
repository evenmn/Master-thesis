#include "system.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/hydrogenorbital.h"
#include "WaveFunctions/padejastrow.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/atomicnucleus.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/randomnormal.h"
#include "InitialWeights/initialweights.h"
#include "InitialWeights/ones.h"
#include "Optimization/optimization.h"
#include "Optimization/gradientdescent.h"

using namespace std;

int main() {
    int numberOfDimensions  = 2;
    int numberOfParticles   = 2;
    int numberOfSteps       = int(1e6);
    double omega            = 1.0;          // Oscillator frequency.
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.0;          // Amount of the total steps used
    bool interaction        = true;

    System* system = new System();
    std::vector<class WaveFunction*> WaveFunctionElements;
    WaveFunctionElements.push_back(new class SimpleGaussian(system, 0));
    WaveFunctionElements.push_back(new class PadeJastrow(system, 1));

    system->setWaveFunction             (WaveFunctionElements);
    system->setInitialState             (new RandomNormal(system, numberOfDimensions, numberOfParticles));
    system->setInitialWeights           (new Ones(system, WaveFunctionElements.size()));
    system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setOptimizer                (new GradientDescent(system));
    system->setEquilibrationFraction    (equilibration);
    system->setInteraction              (interaction);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps);
    return 0;
}
