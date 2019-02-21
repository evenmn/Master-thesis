#include "system.h"

#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/gaussian.h"
#include "WaveFunctions/mlgaussian.h"
#include "WaveFunctions/hydrogenorbital.h"
#include "WaveFunctions/padejastrow.h"
#include "WaveFunctions/nqsjastrow.h"
#include "WaveFunctions/nqsjastrowreal.h"
#include "WaveFunctions/partlyrestricted.h"
#include "WaveFunctions/slaterdeterminant.h"

#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/atomicnucleus.h"

#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/randomnormal.h"

#include "InitialWeights/initialweights.h"
#include "InitialWeights/constant.h"
#include "InitialWeights/randomize.h"

#include "Metropolis/metropolis.h"
#include "Metropolis/bruteforce.h"
#include "Metropolis/importancesampling.h"

#include "Optimization/optimization.h"
#include "Optimization/gradientdescent.h"
#include "Optimization/barzilaiborwein.h"
#include "Optimization/gradientdescentmomentum.h"

using namespace std;

int main() {
    int     numberOfDimensions  = 3;
    int     numberOfParticles   = 1;
    int     numberOfHiddenNodes = 2;
    int     numberOfSteps       = int(1e6);
    int     numberOfIterations  = 50;
    double  eta                 = 0.5;         // Learning rate
    double  omega               = 1.0;          // Oscillator frequency
    double  sigma               = 1.0;          // Width of probability distribution
    double  stepLength          = 0.1;          // Metropolis step length
    double  equilibration       = 0.1;          // Amount of the total steps used
    bool    interaction         = false;
    int     maxNumberOfParametersPerElement = numberOfParticles*numberOfDimensions*numberOfParticles*numberOfDimensions;

    System* system = new System();
    system->setEquilibrationFraction    (equilibration);
    system->setInteraction              (interaction);
    system->setStepLength               (stepLength);
    system->setFrequency                (omega);
    system->setWidth                    (sigma);
    system->setLearningRate             (eta);
    system->setNumberOfParticles        (numberOfParticles);
    system->setNumberOfDimensions       (numberOfDimensions);
    system->setNumberOfHiddenNodes      (numberOfHiddenNodes);
    system->setMaxNumberOfParametersPerElement (maxNumberOfParametersPerElement);
    system->setNumberOfOrbitals         ();
    system->setNumberOfFreeDimensions   ();

    std::vector<class WaveFunction*> WaveFunctionElements;
    WaveFunctionElements.push_back      (new class HydrogenOrbital      (system, 0));
    //WaveFunctionElements.push_back      (new class Gaussian             (system, 0));
    //WaveFunctionElements.push_back      (new class MLGaussian           (system, 0));
    //WaveFunctionElements.push_back      (new class NQSJastrow           (system, 1));
    //WaveFunctionElements.push_back      (new class PartlyRestricted     (system, 1));
    //WaveFunctionElements.push_back      (new class SlaterDeterminant    (system, 1));
    //WaveFunctionElements.push_back      (new class NQSJastrowReal       (system, 1));
    //WaveFunctionElements.push_back      (new class PadeJastrow          (system, 1));

    system->setNumberOfWaveFunctionElements(int(WaveFunctionElements.size()));
    system->setWaveFunction             (WaveFunctionElements);
    system->setInitialWeights           (new Constant(system, 1));
    system->setInitialState             (new RandomNormal(system));
    //system->setHamiltonian              (new HarmonicOscillator(system));
    system->setHamiltonian              (new AtomicNucleus(system, 2));
    system->setMetropolis               (new ImportanceSampling(system));
    system->setOptimization             (new GradientDescent(system));
    system->setGradients                ();
    system->runMetropolisSteps          (numberOfSteps, numberOfIterations);
    return 0;
}
