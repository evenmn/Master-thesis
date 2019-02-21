#include "hydrogenorbital.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

HydrogenOrbital::HydrogenOrbital(System* system, int elementNumber) :
        WaveFunction(system) {
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_elementNumber = elementNumber;
}

double HydrogenOrbital::calculateRadialVectorElement(Eigen::VectorXd positions, int par) {

    double sqrtElementWise = 0;
    for(int d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += positions(par*m_numberOfDimensions + d) * positions(par*m_numberOfDimensions + d);
    }
    return sqrt(sqrtElementWise);
}

Eigen::VectorXd HydrogenOrbital::calculateRadialVector(Eigen::VectorXd positions) {
    Eigen::VectorXd radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        radialVector(i) = calculateRadialVectorElement(positions, i);
    }
    return radialVector;
}

void HydrogenOrbital::initializeArrays(Eigen::VectorXd positions) {
    m_positions       = positions;
    m_radialVector    = calculateRadialVector(positions);
}

void HydrogenOrbital::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions         = m_positions;
    m_positions            = positions;

    m_oldRadialVector      = m_radialVector;
    int particle = int(pRand/m_numberOfDimensions);
    m_radialVector(particle)  = calculateRadialVectorElement(positions, particle);
}

void HydrogenOrbital::resetArrays() {
    m_positions       = m_oldPositions;
    m_radialVector    = m_oldRadialVector;
}

void HydrogenOrbital::updateParameters(Eigen::MatrixXd parameters) {
    m_alpha = parameters(m_elementNumber, 0);
}

double HydrogenOrbital::evaluate() {
    return exp(- m_alpha * m_numberOfParticles * m_radialVector.sum());
}

double HydrogenOrbital::evaluateSqrd() {
    return exp(- 2 * m_alpha * m_numberOfParticles * m_radialVector.sum());
}

double HydrogenOrbital::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    return - m_alpha * m_numberOfParticles;
}

double HydrogenOrbital::computeSecondDerivative() {;
    return 0;
}

Eigen::VectorXd HydrogenOrbital::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    gradients(0) = 0.5 * m_numberOfParticles;
    return gradients;
}

Eigen::VectorXd HydrogenOrbital::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    return gradients;
}
