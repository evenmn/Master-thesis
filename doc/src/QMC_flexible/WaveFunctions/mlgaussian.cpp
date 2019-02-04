#include "mlgaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

MLGaussian::MLGaussian(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber      = elementNumber;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_omega              = m_system->getFrequency();
    m_sigmaSqrd         = m_system->getWidth() * m_system->getWidth();
}

double MLGaussian::evaluate(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    return exp(-double(m_omega * ((particles - (m_parameters.row(m_elementNumber)).head(m_numberOfParticles)).cwiseAbs2()).sum())/(2 * m_sigmaSqrd));
}

double MLGaussian::evaluateSqrd(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    return exp(-double(m_omega * ((particles - (m_parameters.row(m_elementNumber)).head(m_numberOfParticles)).cwiseAbs2()).sum())/m_sigmaSqrd);
}

double MLGaussian::computeFirstDerivative(int k) {
    m_parameters         = m_system->getWeights();
    m_particles          = m_system->getParticles();
    return - m_omega * (m_particles(0, k) - m_parameters(m_elementNumber, k))/m_sigmaSqrd;
}

double MLGaussian::computeSecondDerivative() {;
    return -m_omega * m_numberOfParticles/m_sigmaSqrd;
}

void MLGaussian::computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k) {
    m_radialVector       = m_system->getRadialVector();
    gradients(0) = 0.5 * m_omega * m_radialVector(k);
}

void MLGaussian::computeSecondEnergyDerivative(Eigen::VectorXd &gradients) {
    gradients(0) = 0.5 * m_omega *  m_numberOfParticles * m_numberOfDimensions;
}
