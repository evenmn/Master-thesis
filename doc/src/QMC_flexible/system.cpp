#include "system.h"
#include <cassert>
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Optimization/optimization.h"
#include "Math/random2.h"
#include <iostream>

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    Random2 rand;

    int pRand = rand.nextInt(m_numberOfParticles);
    int dRand = rand.nextInt(m_numberOfDimensions);

    Eigen::MatrixXd newPositions = m_particles;
    newPositions(pRand, dRand) = m_particles(pRand, dRand) + 2*(rand.nextDouble() - 0.5) * m_stepLength;

    //std::cout << m_particles << std::endl;
    //std::cout << " " << std::endl;

    double psiOld = evaluateWaveFunction(m_particles);
    double psiNew = evaluateWaveFunction(newPositions);

    double w = psiNew/psiOld;
    double r = rand.nextDouble();

    //std::cout << w * w << std::endl;
    //std::cout << r << std::endl;
    //std::cout << " " << std::endl;

    if(w * w > r) {
        m_particles(pRand, dRand) = newPositions(pRand, dRand);
        return true;
    }

    return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    int iterations = 1;

    for (int iter = 0; iter < iterations; iter++) {
        for (int i=0; i < numberOfMetropolisSteps; i++) {
            bool acceptedStep = metropolisStep();

            /* Here you should sample the energy (and maybe other things using
             * the m_sampler instance of the Sampler class. Make sure, though,
             * to only begin sampling after you have let the system equilibrate
             * for a while. You may handle this using the fraction of steps which
             * are equilibration steps; m_equilibrationFraction.
             */

            if(double(i)/m_numberOfMetropolisSteps >= m_equilibrationFraction) {
                m_sampler->sample(acceptedStep, i);
            }
        }
        m_sampler->computeAverages();
        m_sampler->printOutputToTerminal();
    }
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(std::vector<class WaveFunction *> waveFunctionVector) {
    m_waveFunctionVector = waveFunctionVector;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setOptimizer(Optimization* optimization) {
    m_optimizer = optimization;
}

double System::evaluateWaveFunction(const Eigen::MatrixXd &particles) {
    double WF = 0;
    for(unsigned i = 0; i < m_waveFunctionVector.size(); i++) {
        WF += m_waveFunctionVector[i]->evaluate(particles);
    }
    return WF;
}

double System::getKineticEnergy(const Eigen::MatrixXd &particles) {
    double KineticEnergy = 0;
    for(unsigned i = 0; i < m_waveFunctionVector.size(); i++) {
        KineticEnergy += m_waveFunctionVector[i]->computeSecondDerivative(particles);
    }
    for(int k = 0; k < m_numberOfParticles; k++) {
        double NablaLnPsi = 0;
        for(unsigned i = 0; i < m_waveFunctionVector.size(); i++) {
            NablaLnPsi += m_waveFunctionVector[i]->computeFirstDerivative(particles, k);
        }
        KineticEnergy += NablaLnPsi * NablaLnPsi;
    }
    return -0.5 * KineticEnergy;
}
