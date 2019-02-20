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

void MLGaussian::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;
}

void MLGaussian::resetArrays() {
    m_positions = m_oldPositions;
}

void MLGaussian::initializeArrays(Eigen::VectorXd positions) {

}

double MLGaussian::evaluate(Eigen::VectorXd positions) {
    m_parameters         = m_system->getWeights();
    //Eigen::VectorXd a    = (m_parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);

    Eigen::VectorXd Xa = positions + Eigen::VectorXd::Ones(m_numberOfFreeDimensions) - m_parameters.row(0).head(m_numberOfFreeDimensions).transpose();
    return exp(-0.5*double(Xa.transpose() * Xa));
    //return exp(-double(m_omega * ((positions + Eigen::VectorXd::Ones(m_numberOfFreeDimensions) - a).cwiseAbs2()).sum())/(2 * m_sigmaSqrd));
}

double MLGaussian::evaluateSqrd(Eigen::VectorXd positions) {
    m_parameters         = m_system->getWeights();
    //Eigen::VectorXd a    = (m_parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);

    Eigen::VectorXd Xa = positions + Eigen::VectorXd::Ones(m_numberOfFreeDimensions) - m_parameters.row(0).head(m_numberOfFreeDimensions).transpose();
    return exp(-double(Xa.transpose() * Xa)/m_sigmaSqrd);

    //return exp(-double(m_omega * ((positions + Eigen::VectorXd::Ones(m_numberOfFreeDimensions)- a).cwiseAbs2()).sum())/m_sigmaSqrd);
}

double MLGaussian::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    m_parameters         = m_system->getWeights();
    return - m_omega * (positions(k) + 1 - m_parameters(m_elementNumber, k))/m_sigmaSqrd;
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
