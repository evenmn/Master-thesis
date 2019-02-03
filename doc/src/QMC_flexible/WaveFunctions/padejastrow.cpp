#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrow::PadeJastrow(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber      = elementNumber;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

double PadeJastrow::evaluate(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            int l = m_numberOfParticles*i + j + 1;      // Stack Gamma matrix
            PadeJastrowFactor += m_parameters(m_elementNumber, l) * distanceMatrix(i,j)/(1 + m_parameters(m_elementNumber, 0) * distanceMatrix(i,j));
        }
    }
    return exp(PadeJastrowFactor);
}

double PadeJastrow::evaluateSqrd(Eigen::MatrixXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    m_parameters         = m_system->getWeights();
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            int l = m_numberOfParticles*i + j + 1;      // Stack Gamma matrix
            PadeJastrowFactor += m_parameters(m_elementNumber, l) * distanceMatrix(i,j)/(1 + m_parameters(m_elementNumber, 0) * distanceMatrix(i,j));
        }
    }
    return exp(2*PadeJastrowFactor);
}

double PadeJastrow::computeFirstDerivative(int k) {
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();
    double derivative = 0;
    for(int j=0; j<k; j++) {
        double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(k,j));
        int l = m_numberOfParticles*k + j + 1;      // Stack Gamma matrix
        derivative+= m_parameters(m_elementNumber, l) * f * f;
    }
    return derivative;
}

double PadeJastrow::computeSecondDerivative() {
    m_radialVector       = m_system->getRadialVector();
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();
    double derivative = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(i,j));
            int l = m_numberOfParticles*i + j + 1;      // Stack Gamma matrix
            derivative += m_parameters(m_elementNumber, l) * f * f * (1/m_radialVector(i) - m_parameters(m_elementNumber, 0) * f);
        }
    }
    return 2 * derivative;
}

void PadeJastrow::computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k) {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += m_particles(i,j) * m_particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }
    gradients(0) = 0.5 * double(r.transpose()*r);
}

void PadeJastrow::computeSecondEnergyDerivative(Eigen::VectorXd &gradients) {

}
