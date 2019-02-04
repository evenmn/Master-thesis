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
            int l = m_numberOfParticles * i + j + 1;      // Stack Gamma matrix
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
            int l = m_numberOfParticles * i + j + 1;      // Stack Gamma matrix
            PadeJastrowFactor += m_parameters(m_elementNumber, l) * distanceMatrix(i,j)/(1 + m_parameters(m_elementNumber, 0) * distanceMatrix(i,j));
        }
    }
    return exp(2 * PadeJastrowFactor);
}

double PadeJastrow::computeFirstDerivative(int k) {
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();
    double derivative = 0;
    for(int j=0; j<k; j++) {
        double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(k,j));
        int l = m_numberOfParticles * k + j + 1;      // Stack Beta matrix
        derivative += m_parameters(m_elementNumber, l) * f * f;
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
            int l = m_numberOfParticles * i + j + 1;      // Stack Beta matrix
            derivative += m_parameters(m_elementNumber, l) * f * f * (1/m_radialVector(i) - m_parameters(m_elementNumber, 0) * f);
        }
    }
    return 2 * derivative;
}

void PadeJastrow::computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k) {
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();

    //Update gamma
    double derivative = 0;
    for(int j=0; j<k; j++) {
        double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(k,j));
        int l = m_numberOfParticles * k + j + 1;      // Stack Beta matrix
        derivative -= m_parameters(m_elementNumber, l) * m_distanceMatrix(k,j) * f * f * f;
    }
    gradients(0) = - derivative;

    //Update Beta matrix
    for(int j=0; j<k; j++) {
        double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(k,j));
        int l = m_numberOfParticles * k + j + 1;      // Stack Beta matrix
        gradients(l) = - 0.5 * f * f;
    }
}

void PadeJastrow::computeSecondEnergyDerivative(Eigen::VectorXd &gradients) {
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_radialVector       = m_system->getRadialVector();
    m_parameters         = m_system->getWeights();

    //Update gamma
    double derivative = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(i,j));
            int l = m_numberOfParticles * i + j + 1;      // Stack Beta matrix
            derivative -= m_parameters(m_elementNumber, l) * m_distanceMatrix(i,j) * f * f * f * (1 / m_distanceMatrix(i,j) + 2 / m_radialVector(i) -3 * m_parameters(m_elementNumber, 0) * f);
        }
    }
    gradients(0) = - derivative;

    //Update Beta matrix
    for(int i=0; i < m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_parameters(m_elementNumber, 0) * m_distanceMatrix(i,j));
            int l = m_numberOfParticles * i + j + 1;      // Stack Beta matrix
            gradients(l) = - f * f * (1/m_radialVector(i) - m_parameters(m_elementNumber, 0) * f);
        }
    }
}
