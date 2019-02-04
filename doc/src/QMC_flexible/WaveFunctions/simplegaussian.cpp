#include "simplegaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

SimpleGaussian::SimpleGaussian(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber      = elementNumber;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_omega              = m_system->getFrequency();
}

double SimpleGaussian::evaluate(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    return exp(-0.5 * m_omega * m_parameters(m_elementNumber, 0) * (particles.cwiseAbs2()).sum());
}

double SimpleGaussian::evaluateSqrd(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters        = m_system->getWeights();
    return exp(- m_omega * m_parameters(m_elementNumber, 0) * (particles.cwiseAbs2()).sum());
}

double SimpleGaussian::computeFirstDerivative(int k) {
    m_parameters         = m_system->getWeights();
    m_radialVector       = m_system->getRadialVector();
    return - m_omega * m_parameters(m_elementNumber, 0) * m_radialVector(k);
}

double SimpleGaussian::computeSecondDerivative() {;
    m_parameters         = m_system->getWeights();
    return - m_omega * m_parameters(m_elementNumber, 0) * m_numberOfParticles * m_numberOfDimensions;
}

void SimpleGaussian::computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k) {
    m_radialVector       = m_system->getRadialVector();
    gradients(0) = 0.5 * m_omega * m_radialVector(k);
}

void SimpleGaussian::computeSecondEnergyDerivative(Eigen::VectorXd &gradients) {
    gradients(0) = 0.5 * m_omega *  m_numberOfParticles * m_numberOfDimensions;
}
