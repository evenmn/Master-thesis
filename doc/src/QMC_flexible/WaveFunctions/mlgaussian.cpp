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

void MLGaussian::updateParameters(Eigen::MatrixXd parameters) {
    m_a = (parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
}

void MLGaussian::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_Xa        = positions + Eigen::VectorXd::Ones(m_numberOfFreeDimensions) - m_a;
}

void MLGaussian::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions  = m_positions;
    m_positions     = positions;
    m_oldXa         = m_Xa;
    m_Xa            = positions + Eigen::VectorXd::Ones(m_numberOfFreeDimensions) - m_a;
}

void MLGaussian::resetArrays() {
    m_positions = m_oldPositions;
    m_Xa        = m_oldXa;
}

double MLGaussian::evaluate() {
    return exp(-0.5 * double(m_Xa.transpose() * m_Xa)/m_sigmaSqrd);
}

double MLGaussian::evaluateSqrd() {
    return exp(-double(m_Xa.transpose() * m_Xa)/m_sigmaSqrd);
}

double MLGaussian::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    return - m_omega * (positions(k) - m_a(k))/m_sigmaSqrd;
}

double MLGaussian::computeSecondDerivative() {;
    return - m_omega * m_numberOfFreeDimensions/m_sigmaSqrd;
}

Eigen::VectorXd MLGaussian::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(k) = - 0.5 * m_omega / m_sigmaSqrd;
    return gradients;
}

Eigen::VectorXd MLGaussian::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
