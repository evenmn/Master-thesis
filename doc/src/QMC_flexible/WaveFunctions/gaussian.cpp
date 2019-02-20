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

void Gaussian::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
}

void Gaussian::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;
}

void Gaussian::resetArrays() {
    m_positions = m_oldPositions;
}

void Gaussian::updateParameters(Eigen::MatrixXd parameters) {
    m_alpha = parameters(m_elementNumber,0);
}

double Gaussian::evaluate() {
    return exp(-0.5 * m_omega * m_alpha * (m_positions.cwiseAbs2()).sum());
}

double Gaussian::evaluateSqrd() {
    return exp(- m_omega * m_alpha * (m_positions.cwiseAbs2()).sum());
}

double Gaussian::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    return - m_omega * m_alpha * positions(k);
}

double Gaussian::computeSecondDerivative() {;
    return - m_omega * m_alpha * m_numberOfFreeDimensions;
}

Eigen::VectorXd Gaussian::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = 0.5 * m_omega * m_positions(k);
    return gradients;
}

Eigen::VectorXd Gaussian::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = 0.5 * m_omega *  m_numberOfFreeDimensions;
    return gradients;
}
