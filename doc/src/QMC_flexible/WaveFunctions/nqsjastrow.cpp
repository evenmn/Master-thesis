#include "nqsjastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

NQSJastrow::NQSJastrow(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber                     = elementNumber;
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
}

double NQSJastrow::evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    double PadeJastrowFactor = 0;
    return exp(PadeJastrowFactor);
}

double NQSJastrow::evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    double PadeJastrowFactor = 0;
    return exp(2 * PadeJastrowFactor);
}

double NQSJastrow::computeFirstDerivative(const Eigen::VectorXd particles, int k) {
    m_distanceMatrix = m_system->getDistanceMatrix();
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    return derivative;
}

double NQSJastrow::computeSecondDerivative() {
    m_radialVector       = m_system->getRadialVector();
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();

    double derivative = 0;
    return derivative;
}

Eigen::VectorXd NQSJastrow::computeFirstEnergyDerivative(int k) {
    m_particles = m_system->getParticles();

    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    return gradients;
}

Eigen::VectorXd NQSJastrow::computeSecondEnergyDerivative() {
    m_distanceMatrix     = m_system->getDistanceMatrix();

    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    return gradients;
}
