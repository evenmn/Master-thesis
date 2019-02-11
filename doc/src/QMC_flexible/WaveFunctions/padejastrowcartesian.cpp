#include "padejastrowcartesian.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrowCartesian::PadeJastrowCartesian(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber                     = elementNumber;
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
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
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        double F = f(k_p,j_p);
        derivative += beta(k_p,j_p) * F * F * g(k,j, k_p, j_p);
    }
    return derivative;
}

double PadeJastrowCartesian::computeSecondDerivative() {
    m_radialVector       = m_system->getRadialVector();
    m_distanceMatrix     = m_system->getDistanceMatrix();
    m_parameters         = m_system->getWeights();

    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double F = f(i_p,j_p);
            double G = g(i,j,i_p,j_p);
            derivative += beta(i_p,j_p) * F * F * (1 - G * G * (1 + 3 * gamma() * m_distanceMatrix(i_p,j_p)) * F) / m_distanceMatrix(i_p,j_p);
        }
    }
    return derivative;
}

Eigen::VectorXd PadeJastrowCartesian::computeFirstEnergyDerivative(int k) {
    m_particles = m_system->getParticles();

    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    //Update beta matrix
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        int l = m_numberOfDimensions * k_p + j_p + 1;
        double F = f(k_p,j_p);
        gradients(l) = F * F * g(k,j,k_p,j_p);
    }

    //Update gamma
    double derivative = 0;
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        double F = f(k_p,j_p);
        derivative -= beta(k_p,j_p) * F * F * F * (m_particles(k) - m_particles(j));
    }
    gradients(0) = 2 * derivative;
    return gradients;
}

Eigen::VectorXd PadeJastrowCartesian::computeSecondEnergyDerivative() {
    m_distanceMatrix     = m_system->getDistanceMatrix();

    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    //Update Beta matrix
    for(int i=0; i < m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double F = f(i_p,j_p);
            double G = g(i,j,i_p,j_p);
            int l = m_numberOfParticles * i_p + j_p + 1;      // Stack Beta matrix
            gradients(l) = 0.5 * F * F * (G * G * (1 + 3 * gamma() * m_distanceMatrix(i_p,j_p)) * F - 1) / m_distanceMatrix(i_p,j_p);
        }
    }

    //Update gamma
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double F = f(i_p,j_p);
            double G = g(i,j,i_p,j_p);
            derivative -= beta(i_p,j_p) * F * F * F * (3 * G * G * F * gamma() * m_distanceMatrix(i_p,j_p) - 1);
        }
    }
    gradients(0) = derivative;
    return gradients;
}
