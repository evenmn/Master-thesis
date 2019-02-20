#include "system.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "InitialWeights/initialweights.h"
#include "Metropolis/metropolis.h"
#include "Optimization/optimization.h"
#include "Math/random2.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>
#include <string>

void System::runMetropolisSteps(int numberOfMetropolisSteps, int numberOfIterations) {
    m_positions                 = m_initialState->getParticles();
    m_parameters                = m_initialWeights->getWeights();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    std::string path = "../../data/";                       //Path to data folder
    std::string energy_filename = generate_filename("Energy3", ".dat");
    std::ofstream energy;
    energy.open(path + energy_filename);                    //Open energy file based on parameters
    for (int iter = 0; iter < numberOfIterations; iter++) {
        clock_t start_time = clock();
        for (int i=0; i < numberOfMetropolisSteps; i++) {
            bool acceptedStep = m_metropolis->acceptMove();
            m_positions       = m_metropolis->updatePositions();
            if(double(i)/m_numberOfMetropolisSteps >= m_equilibrationFraction) {
                m_sampler->sample(acceptedStep, i);
            }
        }
        clock_t end_time = clock();
        m_sampler->computeAverages();
        m_sampler->printOutputToTerminal(iter, double(end_time - start_time)/CLOCKS_PER_SEC);
        energy << m_sampler->getEnergy() << "\n";
        m_parameters -= m_optimization->updateParameters();
        updateAllParameters(m_parameters);
    }
    if(energy.is_open())  energy.close();
}

int factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

double binomial(int n, int p) {
    //Binomial coefficients, equal to magic numbers
    return factorial(n+p)/(factorial(n)*factorial(p));
}

int orbitals(int numberOfParticles, int numberOfDimensions) {

    int counter = 0;
    while(true) {
        double orb = 2*binomial(counter, numberOfDimensions);
        if(int(orb) == numberOfParticles) {
            return counter+1;
            break;
        }
        else if(orb > numberOfParticles) {
            std::cout << "This program supports closed-shells only. Please choose a P such that the orbital is full" << std::endl;
            exit(0);
        }
        counter += 1;
    }
}

void System::setNumberOfOrbitals() {
    m_numberOfOrbitals = orbitals(m_numberOfParticles, m_numberOfDimensions);
}

void System::setNumberOfParticles(int numberOfParticles) {
    assert(numberOfParticles > 0);
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    assert(numberOfDimensions > 0);
    m_numberOfDimensions = numberOfDimensions;
}

void System::setNumberOfFreeDimensions() {
    m_numberOfFreeDimensions = m_numberOfParticles * m_numberOfDimensions;
}

void System::setNumberOfHiddenNodes(int numberOfHiddenNodes) {
    assert(numberOfHiddenNodes > 0);
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

void System::setNumberOfWaveFunctionElements(int numberOfWaveFunctionElements) {
    m_numberOfWaveFunctionElements = numberOfWaveFunctionElements;
}

void System::setMaxNumberOfParametersPerElement(int maxNumberOfParametersPerElement) {
    m_maxNumberOfParametersPerElement = maxNumberOfParametersPerElement;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setInteraction(bool interaction) {
    m_interaction = interaction;
}

void System::setFrequency(double omega) {
    assert(omega > 0);
    m_omega = omega;
}

void System::setLearningRate(double eta) {
    assert(eta > 0);
    m_eta = eta;
}

void System::setWidth(double sigma) {
    assert(sigma > 0);
    m_sigma = sigma;
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

void System::setInitialWeights(InitialWeights* initialWeights) {
    m_initialWeights = initialWeights;
}

void System::setMetropolis(Metropolis* metropolis) {
    m_metropolis = metropolis;
}

void System::setOptimization(Optimization* optimization) {
    m_optimization = optimization;
}

void System::setGradients() {
    int maxNumberOfParametersPerElement = m_numberOfParticles*m_numberOfParticles + m_numberOfParticles;
    m_gradients = Eigen::MatrixXd(m_waveFunctionVector.size(), maxNumberOfParametersPerElement);
}

void System::updateAllArrays(Eigen::VectorXd particles, int pRand) {
    for(auto& i : m_waveFunctionVector) {
        i->updateArrays(particles, pRand);
    }
}

void System::resetAllArrays() {
    for(auto& i : m_waveFunctionVector) {
        i->resetArrays();
    }
}

void System::updateAllParameters(Eigen::MatrixXd parameters) {
    for(auto& i : m_waveFunctionVector) {
        i->updateParameters(parameters);
    }
}

double System::evaluateWaveFunction() {
    double WF = 1;
    for(auto& i : m_waveFunctionVector) {
        WF *= i->evaluate();
    }
    return WF;
}

double System::evaluateWaveFunctionSqrd() {
    double WF = 1;
    for(auto& i : m_waveFunctionVector) {
        WF *= i->evaluateSqrd();
    }
    return WF;
}

double System::getKineticEnergy() {
    double KineticEnergy = 0;
    for(auto& i : m_waveFunctionVector) {
        KineticEnergy += i->computeSecondDerivative();
    }
    for(int k = 0; k < m_numberOfFreeDimensions; k++) {
        double NablaLnPsi = 0;
        for(auto& i : m_waveFunctionVector) {
            NablaLnPsi += i->computeFirstDerivative(m_positions, k);
        }
        KineticEnergy += NablaLnPsi * NablaLnPsi;
    }
    return -0.5 * KineticEnergy;
}

std::string System::generate_filename(std::string name, std::string extension) {
    return name + extension;
}
