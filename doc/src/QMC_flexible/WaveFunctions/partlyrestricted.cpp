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
    m_omega                             = m_system->getFrequency();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd = sigma*sigma;
}

void PartlyRestricted::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;
}

void PartlyRestricted::resetArrays() {
    m_positions = m_oldPositions;
}

void PartlyRestricted::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
}

double PartlyRestricted::evaluate(Eigen::VectorXd positions) {
    m_parameters         = m_system->getWeights();

    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            double c = m_parameters(m_elementNumber, j*m_numberOfFreeDimensions + i);
            Sum += positions(i) * positions(j) * c;
        }
    }
    return exp(-0.5 * Sum  / (m_sigmaSqrd * m_sigmaSqrd));
}

double PartlyRestricted::evaluateSqrd(Eigen::VectorXd positions) {
    m_parameters         = m_system->getWeights();
    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            double c = m_parameters(m_elementNumber, j*m_numberOfFreeDimensions + i);
            Sum += positions(i) * positions(j) * c;
        }
    }
    return exp(- Sum  / (m_sigmaSqrd * m_sigmaSqrd));
}

double PartlyRestricted::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    m_parameters         = m_system->getWeights();
    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        double c = m_parameters(m_elementNumber, k*m_numberOfFreeDimensions + i);
        Sum += positions(i) * c;
    }
    return - Sum / (m_sigmaSqrd + m_sigmaSqrd);
}

double PartlyRestricted::computeSecondDerivative() {
    m_parameters         = m_system->getWeights();
    double Sum = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        Sum += m_parameters(m_elementNumber, i*m_numberOfFreeDimensions + i);
    }
    return - Sum / (m_sigmaSqrd + m_sigmaSqrd);
}

Eigen::VectorXd PartlyRestricted::computeFirstEnergyDerivative(int k) {
    m_positions = m_system->getPositions();
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        gradients(k * m_numberOfFreeDimensions + i) = 0.5 * m_positions(i) / (m_sigmaSqrd + m_sigmaSqrd);
    }
    return gradients;
}

Eigen::VectorXd PartlyRestricted::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        gradients(i * m_numberOfFreeDimensions + i) = 0.5 / (m_sigmaSqrd + m_sigmaSqrd);
    }
    return gradients;
}
