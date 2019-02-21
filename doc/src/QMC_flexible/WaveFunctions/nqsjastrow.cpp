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

void NQSJastrow::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldV = m_v;
    m_oldN = m_n;
    m_oldP = m_p;
    m_v = m_b + m_W.transpose() * positions;

    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        m_n(i) = 1/(1 + exp(-m_v(i)));
        m_p(i) = 1/(1 + exp(+m_v(i)));
    }
}

void NQSJastrow::resetArrays() {
    m_positions = m_oldPositions;
    m_v         = m_oldV;
    m_n         = m_oldN;
    m_p         = m_oldP;
}

void NQSJastrow::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_v = m_b + m_W.transpose() * positions;

    for(int i=0; i<m_numberOfHiddenNodes; i++) {
        m_n(i) = 1/(1 + exp(-m_v(i)));
        m_p(i) = 1/(1 + exp(+m_v(i)));
    }
}

void NQSJastrow::updateParameters(Eigen::MatrixXd parameters) {
    Eigen::VectorXd XXX = parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(XXX.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    m_W = W;

    m_b = parameters.row(m_elementNumber).head(m_numberOfHiddenNodes);
}

int fromWToParameterIndex(int i, int j, int numberOfFreeDimensions) {
    return j*numberOfFreeDimensions + i;
}

double NQSJastrow::evaluate() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= 1 + exp(m_v(j));
    }
    return Prod;
}

double NQSJastrow::evaluateSqrd() {
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        Prod *= 1 + exp(m_v(j));
    }
    return Prod * Prod;
}

double NQSJastrow::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    double Sum1 = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        double Sum2 = 0;
        for(int i=0; i<m_numberOfFreeDimensions; i++) {
            Sum2 += m_W(i,j) * positions(i) / m_sigmaSqrd;
        }
        Sum1 += m_W(k,j) / (m_sigmaSqrd * (1 + exp(-m_b(j) - Sum2)));
    }
    return Sum1;
}

double NQSJastrow::computeSecondDerivative() {
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_W(k,j) * m_W(k,j) * m_n(j) * m_p(j);
        }
    }
    return Sum / (m_sigmaSqrd * m_sigmaSqrd);
}

Eigen::VectorXd NQSJastrow::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_W(k,j) * m_n(j) * m_p(j) / m_sigmaSqrd;
        }
    }
    gradients(0) = -0.5 * Sum;


    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int n = fromWToParameterIndex(l, m, m_numberOfFreeDimensions);
            if(l == k) {
                gradients(n + m_numberOfHiddenNodes) = -0.5 * m_n(m) * (m_sigmaSqrd + m_W(k,m) * m_p(m) * m_positions(k));
            }
            else {
                gradients(n + m_numberOfHiddenNodes) = -0.5 * m_W(k, m) * m_n(m) * m_p(m) * m_positions(l);
            }
        }
    }
    return gradients / (m_sigmaSqrd * m_sigmaSqrd);
}

Eigen::VectorXd NQSJastrow::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    double Sum = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            Sum += m_W(k,j) * m_W(k,j) * m_n(j) * m_p(j) * (m_p(j) - m_n(j));
        }
    }
    gradients(0) = -0.5 * Sum / (m_sigmaSqrd * m_sigmaSqrd);

    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int k = fromWToParameterIndex(l, m, m_numberOfFreeDimensions);
            gradients(k + m_numberOfHiddenNodes) = m_n(m) * m_p(m) * (2 * m_W(l,m) + m_W.cwiseAbs2().colwise().sum()(m) * m_positions(l) * (m_p(m) - m_n(m)));
        }
    }

    return gradients / (m_sigmaSqrd * m_sigmaSqrd);
}
