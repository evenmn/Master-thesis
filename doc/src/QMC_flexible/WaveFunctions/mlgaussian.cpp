#include "mlgaussian.h"
#include <cassert>
#include <iostream>
#include "../system.h"

MLGaussian::MLGaussian(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber                     = elementNumber;
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd = sigma*sigma;
}

double MLGaussian::evaluate(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    Eigen::VectorXd a    = (m_parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
    return exp(-double(m_omega * ((positions - a).cwiseAbs2()).sum())/(2 * m_sigmaSqrd));
}

double MLGaussian::evaluateSqrd(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    Eigen::VectorXd a    = (m_parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
    return exp(-double(m_omega * ((positions - a).cwiseAbs2()).sum())/m_sigmaSqrd);
}

double MLGaussian::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    m_parameters         = m_system->getWeights();
    return - m_omega * (positions(k) - m_parameters(m_elementNumber, k))/m_sigmaSqrd;
}

double MLGaussian::computeSecondDerivative() {;
    return - m_omega * m_numberOfFreeDimensions/m_sigmaSqrd;
}

Eigen::VectorXd MLGaussian::computeFirstEnergyDerivative(int k) {
    return - (0.5 * m_omega / m_sigmaSqrd) * Eigen::VectorXd::Ones(m_maxNumberOfParametersPerElement);
}

Eigen::VectorXd MLGaussian::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
