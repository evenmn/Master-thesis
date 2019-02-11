#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "sampler.h"
#include "system.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
//#include "Optimization/optimization.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber         = 0;
    m_acceptenceRatio    = 0;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfElements   = m_system->getNumberOfWaveFunctionElements();
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep, int stepNumber) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (stepNumber == int(m_numberOfMetropolisSteps * m_system->getEquilibrationFraction())) {
        m_acceptenceRatio = 0;
        m_cumulativeEnergy = 0;
        int maxNumberOfParametersPerElement = m_numberOfParticles * m_numberOfParticles + m_numberOfParticles;
        m_dE = Eigen::MatrixXd::Zero(m_numberOfElements, maxNumberOfParametersPerElement);
        m_dEE = Eigen::MatrixXd::Zero(m_numberOfElements, maxNumberOfParametersPerElement);
        m_SqrdE = 0;
    }
    m_stepNumber = stepNumber;

    int maxNumberOfParametersPerElement = m_numberOfParticles * m_numberOfParticles + m_numberOfParticles;
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfElements, maxNumberOfParametersPerElement);

    double EL = m_system->getHamiltonian()->computeLocalEnergy();
    m_system->updateParameters(gradients);

    m_cumulativeEnergy  += EL;
    m_dE += gradients;
    m_dEE += gradients * EL;
    m_SqrdE += EL * EL;

    if(acceptedStep) { m_acceptenceRatio += 1; }
    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = 1; //m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    //std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << m_numberOfParticles << endl;
    cout << " Number of dimensions : " << m_numberOfDimensions << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    //for (int i=0; i < p; i++) {
    //    cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    //}
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy           : " << m_energy << endl;
    cout << " Acceptence Ratio : " << double(m_acceptenceRatio)/m_numberOfMetropolisSteps << endl;
    cout << " Variance         : " << 2*(m_SqrdE/m_numberOfMetropolisSteps - m_energy * m_energy) << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     */
    m_energy = m_cumulativeEnergy / ((1 - m_system->getEquilibrationFraction()) * m_system->getNumberOfMetropolisSteps());
}

Eigen::MatrixXd Sampler::getEnergyGradient() {
    double eta = m_system->getLearningRate();
    return 2 * eta * (m_dEE - m_energy * m_dE)/m_system->getNumberOfMetropolisSteps();
}
