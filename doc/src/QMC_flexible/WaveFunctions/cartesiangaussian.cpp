#include "cartesiangaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

CartesianGaussian::CartesianGaussian(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber      = elementNumber;
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_omega              = m_system->getFrequency();
    //m_alpha              = (m_system->getWeights())(m_elementNumber,0);
}

double CartesianGaussian::evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return exp(-0.5 * m_omega * m_alpha * (particles.cwiseAbs2()).sum());
}

double CartesianGaussian::evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return exp(- m_omega * m_alpha * (particles.cwiseAbs2()).sum());
}

double CartesianGaussian::computeFirstDerivative(int k) {
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    m_particles          = m_system->getParticles();
    return - m_omega * m_alpha * m_particles(k);
}

double CartesianGaussian::computeSecondDerivative() {;
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return - m_omega * m_alpha * m_numberOfFreeDimensions;
}

Eigen::VectorXd CartesianGaussian::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    m_particles          = m_system->getParticles();
    gradients(0) = 0.5 * m_omega * m_particles(k);
    return gradients;
}

Eigen::VectorXd CartesianGaussian::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = 0.5 * m_omega *  m_numberOfFreeDimensions;
    return gradients;
}
