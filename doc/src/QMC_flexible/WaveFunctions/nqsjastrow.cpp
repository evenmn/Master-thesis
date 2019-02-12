#include "nqsjastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

NQSJastrow::NQSJastrow(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber                     = elementNumber;
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
}

Eigen::MatrixXd NQSJastrow::W() {
    m_parameters = m_system->getWeights();
    Eigen::VectorXd XXX = m_parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(XXX.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    return W;
}

Eigen::VectorXd NQSJastrow::b() {
    m_parameters = m_system->getWeights();
    return (m_parameters.row(m_elementNumber)).head(m_numberOfHiddenNodes);
}

Eigen::VectorXd NQSJastrow::v(Eigen::VectorXd positions) {
    return b() + W().transpose() * positions;
}

Eigen::VectorXd NQSJastrow::f(Eigen::VectorXd positions) {
    Eigen::VectorXd V = v(positions);
    Eigen::VectorXd F = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        F(i) = exp(V(i));
    }
    return F;
}

Eigen::VectorXd NQSJastrow::g(Eigen::VectorXd positions) {
    Eigen::VectorXd V = v(positions);
    Eigen::VectorXd G = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);
    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        G(i) = 1/(1 + exp(V(i)));
    }
    return G;
}

double NQSJastrow::evaluate(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    Eigen::VectorXd G = g(positions);
    return (G.cwiseInverse()).prod();
}

double NQSJastrow::evaluateSqrd(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    Eigen::VectorXd G = g(positions);
    return ((G.cwiseInverse()).cwiseAbs2()).prod();
}

double NQSJastrow::computeFirstDerivative(const Eigen::VectorXd positions, int k) {
    Eigen::VectorXd F = f(positions);
    Eigen::VectorXd G = g(positions);
    Eigen::MatrixXd w = W();
    return double(w.row(k) * (F.cwiseProduct(G))) / m_sigmaSqrd;
}

double NQSJastrow::computeSecondDerivative() {
    m_positions = m_system->getParticles();
    Eigen::VectorXd F = f(m_positions);
    Eigen::VectorXd G = g(m_positions);
    Eigen::MatrixXd w = W();
    return (w.cwiseAbs2() * F.cwiseProduct(G.cwiseAbs2())).sum();
}

Eigen::VectorXd NQSJastrow::computeFirstEnergyDerivative(int k) {
    m_positions = m_system->getParticles();
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    return gradients;
}

Eigen::VectorXd NQSJastrow::computeSecondEnergyDerivative() {
    m_distanceMatrix     = m_system->getDistanceMatrix();
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    return gradients;
}
