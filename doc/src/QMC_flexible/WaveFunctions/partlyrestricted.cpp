#include "partlyrestricted.h"
#include <cassert>
#include <iostream>
#include "../system.h"

PartlyRestricted::PartlyRestricted(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber                     = elementNumber;
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd2 = sigma*sigma*sigma*sigma;
}

void PartlyRestricted::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldXCx = m_xCx;
    m_xCx = positions.transpose() * m_c * positions;
}

void PartlyRestricted::resetArrays() {
    m_positions = m_oldPositions;
    m_xCx       = m_oldXCx;
}

void PartlyRestricted::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_xCx = positions.transpose() * m_c * positions;
}

void PartlyRestricted::updateParameters(Eigen::MatrixXd parameters) {
    Eigen::Map<Eigen::MatrixXd> c(parameters.row(m_elementNumber).data(), m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    m_c = c;
}

double PartlyRestricted::evaluate() {
    return exp(-0.5 * m_xCx  / m_sigmaSqrd2);
}

double PartlyRestricted::evaluateSqrd() {
    return exp(- m_xCx  / m_sigmaSqrd2);
}

double PartlyRestricted::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    return - double(m_c.row(k) * positions) / m_sigmaSqrd2;
}

double PartlyRestricted::computeSecondDerivative() {
    return - m_c.diagonal().sum() / m_sigmaSqrd2;
}

Eigen::VectorXd PartlyRestricted::computeFirstEnergyDerivative(int k) {
    m_positions = m_system->getPositions();
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        gradients(k * m_numberOfFreeDimensions + i) = 0.5 * m_positions(i) / m_sigmaSqrd2;
    }
    return gradients;
}

Eigen::VectorXd PartlyRestricted::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        gradients(i * m_numberOfFreeDimensions + i) = 0.5 / m_sigmaSqrd2;
    }
    return gradients;
}
