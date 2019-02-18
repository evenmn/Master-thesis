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
    m_radialVector              = m_initialState->getRadialVector();
    m_distanceMatrix            = m_initialState->getDistanceMatrix();
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
            m_radialVector    = m_metropolis->updateRadialVector();
            m_distanceMatrix  = m_metropolis->updateDistanceMatrix();
            if(double(i)/m_numberOfMetropolisSteps >= m_equilibrationFraction) {
                m_sampler->sample(acceptedStep, i);
            }
        }
        clock_t end_time = clock();
        m_sampler->computeAverages();
        m_sampler->printOutputToTerminal(iter, double(end_time - start_time)/CLOCKS_PER_SEC);
        energy << m_sampler->getEnergy() << "\n";
        m_parameters -= m_optimization->updateParameters();
    }
    if(energy.is_open())  energy.close();
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

double System::evaluateWaveFunction(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    double WF = 1;
    for(auto& i : m_waveFunctionVector) {
        WF *= i->evaluate(particles, radialVector, distanceMatrix);
    }
    return WF;
}

double System::evaluateWaveFunctionSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    double WF = 1;
    for(auto& i : m_waveFunctionVector) {
        WF *= i->evaluateSqrd(particles, radialVector, distanceMatrix);
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

Eigen::VectorXd System::calculateRadialVector(Eigen::VectorXd particles) {
    Eigen::VectorXd radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int d=0; d<m_numberOfDimensions; d++) {
            sqrtElementWise += particles(i*m_numberOfDimensions + d) * particles(i*m_numberOfDimensions + d);
        }
        radialVector(i) = sqrt(sqrtElementWise);
    }
    return radialVector;
}

Eigen::MatrixXd System::calculateDistanceMatrix(Eigen::VectorXd particles) {
    Eigen::MatrixXd distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double sqrtElementWise = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                double numb = particles(i*m_numberOfDimensions + d) - particles(j*m_numberOfDimensions + d);
                sqrtElementWise += numb * numb;
            }
            distanceMatrix(i,j) = sqrt(sqrtElementWise);
        }
    }
    return distanceMatrix;
}

double System::calculateDistanceMatrixElement(int i, int j, Eigen::VectorXd particles) {
    // Update element (i,j) in Dist_inv_invance matrix

    double dist = 0;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = particles(m_numberOfDimensions*i+d)-particles(m_numberOfDimensions*j+d);
        dist += diff*diff;
    }
    return dist;
}


void System::calculateDistanceMatrixCross(int par, Eigen::VectorXd particles, Eigen::MatrixXd &distanceMatrix) {
    // Update distance matrix when position of particle "par" is changed

    // Update row
    for(int i=0; i<par; i++) {
        distanceMatrix(par, i) = calculateDistanceMatrixElement(par, i, particles);
    }
    // Update column
    for(int i=par+1; i<m_numberOfParticles; i++) {
        distanceMatrix(i, par) = calculateDistanceMatrixElement(i, par, particles);
    }
}

std::string System::generate_filename(std::string name, std::string extension) {
    return name + extension;
}
