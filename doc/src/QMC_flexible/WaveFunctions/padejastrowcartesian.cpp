#include "padejastrowcartesian.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrowCartesian::PadeJastrowCartesian(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber      = elementNumber;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

double PadeJastrowCartesian::f(int i, int j) {
    m_radialVector = m_system->getRadialVector();
    m_parameters   = m_system->getWeights();
    return 1/(1 + m_parameters(m_elementNumber, 0) * m_radialVector(i,j));
}

double PadeJastrowCartesian::g(int i, int j, int k, int l) {
    m_distanceMatrix    = m_system->getDistanceMatrix();
    m_particles         = m_system->getParticles();
    return (m_particles(i) - m_particles(j))/m_distanceMatrix(k,l);
}

double PadeJastrowCartesian::beta(int i, int j) {
    m_parameters  = m_system->getWeights();
    int k = m_numberOfDimensions * i + j + 1;      // Stack Gamma matrix
    return m_parameters(m_elementNumber, k);
}

double PadeJastrowCartesian::gamma() {
    m_parameters = m_system->getWeights();
    return m_parameters(m_elementNumber, 0);
}

double PadeJastrowCartesian::evaluate(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + gamma() * radialVector(i,j));
            PadeJastrowFactor += beta(i,j) * f * distanceMatrix(i,j);
        }
    }
    return exp(PadeJastrowFactor);
}

double PadeJastrowCartesian::evaluateSqrd(Eigen::VectorXd particles, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + gamma() * radialVector(i,j));
            PadeJastrowFactor += beta(i,j) * f * distanceMatrix(i,j);
        }
    }
    return exp(2 * PadeJastrowFactor);
}

double PadeJastrowCartesian::computeFirstDerivative(int k) {
    int k_particle = int(k/m_numberOfDimensions);

    double derivative = 0;
    for(int j=0; j<k_particle; j++) {
        int j_particle = int(j/m_numberOfDimensions);
        derivative += beta(k_particle,j_particle) * f(k_particle,j_particle) * f(k_particle,j_particle) * g(k,j, k_particle, j_particle);
    }
    return derivative;
}

double PadeJastrowCartesian::computeSecondDerivative() {
    m_radialVector       = m_system->getRadialVector();
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();

    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_particle = int(i/m_numberOfDimensions);
        for(int j=0; j<i_particle; j++) {
            int j_particle = int(j/m_numberOfDimensions);
            derivative += beta(i_particle,j_particle) * f(i_particle,j_particle) * f(i_particle,j_particle) * (1 - g(i, j, i_particle,j_particle) * g(i, j, i_particle,j_particle) * (1 + 3 * gamma() * m_distanceMatrix(i_particle,j_particle)) * f(i_particle,j_particle)) / m_distanceMatrix(i_particle,j_particle);
        }
    }
    return derivative;
}

void PadeJastrowCartesian::computeFirstEnergyDerivative(Eigen::VectorXd &gradients, int k) {
    m_particles = m_system->getParticles();

    int k_particle = int(k/m_numberOfDimensions);

    //Update beta matrix
    for(int j=0; j<k_particle; j++) {
        int j_particle = int(j/m_numberOfDimensions);
        int l = m_numberOfDimensions * k_particle + j_particle + 1;
        gradients(l) = f(k_particle,j_particle) * f(k_particle,j_particle) * g(k, j, k_particle,j_particle);
    }

    //Update gamma
    double derivative = 0;
    for(int j=0; j<k_particle; j++) {
        int j_particle = int(j/m_numberOfDimensions);
        derivative -= beta(k_particle,j_particle) * f(k_particle,j_particle) * f(k_particle,j_particle) * f(k_particle,j_particle) * (m_particles(k) - m_particles(j));
    }
    gradients(0) = 2 * derivative;
}

void PadeJastrowCartesian::computeSecondEnergyDerivative(Eigen::VectorXd &gradients) {
    m_distanceMatrix     = m_system->getDistanceMatrix();

    //Update Beta matrix
    for(int i=0; i < m_numberOfFreeDimensions; i++) {
        int i_particle = int(i/m_numberOfDimensions);
        for(int j=0; j<i; j++) {
            int j_particle = int(j/m_numberOfDimensions);
            int l = m_numberOfParticles * i_particle + j_particle + 1;      // Stack Beta matrix
            gradients(l) = 0.5 * f(i_particle,j_particle) * f(i_particle,j_particle) * (g(i_particle,j_particle, i, j) * g(i, j, i_particle, j_particle) * (1 + 3 * gamma() * m_distanceMatrix(i_particle,j_particle)) * f(i_particle,j_particle) - 1) / m_distanceMatrix(i_particle,j_particle);
        }
    }

    //Update gamma
    double derivative = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        int i_particle = int(i/m_numberOfDimensions);
        for(int j=0; j<i; j++) {
            int j_particle = int(j/m_numberOfDimensions);
            derivative -= beta(i_particle,j_particle) * f(i_particle,j_particle) * f(i_particle,j_particle) * f(i_particle,j_particle) * (3 * g(i,j,i_particle, j_particle) * g(i,j,i_particle, j_particle) * f(i_particle,j_particle) * gamma() * m_distanceMatrix(i_particle,j_particle) - 1);
        }
    }
    gradients(0) = derivative;
}
