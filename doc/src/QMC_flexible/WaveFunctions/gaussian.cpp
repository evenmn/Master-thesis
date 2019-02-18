#include "gaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

Gaussian::Gaussian(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber      = elementNumber;
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_omega              = m_system->getFrequency();
}

double Gaussian::evaluate(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return exp(-0.5 * m_omega * m_alpha * (positions.cwiseAbs2()).sum());
}

double Gaussian::evaluateSqrd(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return exp(- m_omega * m_alpha * (positions.cwiseAbs2()).sum());
}

double Gaussian::computeFirstDerivative(const Eigen::VectorXd positions, int k) {
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return - m_omega * m_alpha * positions(k);
}

double Gaussian::computeSecondDerivative() {;
    m_alpha              = (m_system->getWeights())(m_elementNumber,0);
    return - m_omega * m_alpha * m_numberOfFreeDimensions;
}

Eigen::VectorXd Gaussian::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    m_positions          = m_system->getParticles();
    gradients(0) = 0.5 * m_omega * m_positions(k);
    return gradients;
}

Eigen::VectorXd Gaussian::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = 0.5 * m_omega *  m_numberOfFreeDimensions;
    return gradients;
}
