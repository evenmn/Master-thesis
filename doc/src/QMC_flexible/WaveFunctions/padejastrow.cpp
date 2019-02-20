#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrow::PadeJastrow(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber                     = elementNumber;
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
}

Eigen::MatrixXd PadeJastrow::calculateDistanceMatrix(Eigen::VectorXd particles) {
    Eigen::MatrixXd distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double sqrtElementWise = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                double numb = particles(i*m_numberOfDimensions + d) - particles(j*m_numberOfDimensions + d);
                sqrtElementWise += numb * numb;
            }
            distanceMatrix(i,j) = sqrt(sqrtElementWise);
        }
    }
    return distanceMatrix;
}

double PadeJastrow::calculateDistanceMatrixElement(int i, int j, Eigen::VectorXd particles) {
    // Update element (i,j) in Dist_inv_invance matrix

    double dist = 0;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = particles(m_numberOfDimensions*i+d)-particles(m_numberOfDimensions*j+d);
        dist += diff*diff;
    }
    return dist;
}


Eigen::MatrixXd PadeJastrow::calculateDistanceMatrixCross(int par, Eigen::VectorXd particles) {
    // Update distance matrix when position of particle "par" is changed

    Eigen::MatrixXd distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    // Update row
    for(int i=0; i<par; i++) {
        distanceMatrix(par, i) = calculateDistanceMatrixElement(par, i, particles);
    }
    // Update column
    for(int i=par+1; i<m_numberOfParticles; i++) {
        distanceMatrix(i, par) = calculateDistanceMatrixElement(i, par, particles);
    }
    return distanceMatrix;
}

void PadeJastrow::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_distanceMatrix = calculateDistanceMatrix(positions);
}

void PadeJastrow::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldDistanceMatrix = m_distanceMatrix;
    //m_distanceMatrix = calculateDistanceMatrixCross(int(pRand/m_numberOfDimensions), positions);
    m_distanceMatrix = calculateDistanceMatrix(positions);
}

void PadeJastrow::resetArrays() {
    m_positions = m_oldPositions;
    m_distanceMatrix = m_oldDistanceMatrix;
}


//double PadeJastrowCartesian::f(int i, int j) {
//    m_radialVector = m_system->getRadialVector();
//    m_parameters   = m_system->getWeights();
//    return 1/(1 + m_parameters(m_elementNumber, 0) * m_radialVector(i,j));
//}

Eigen::MatrixXd PadeJastrow::f() {
    m_positions = m_system->getPositions();
    m_distanceMatrix = calculateDistanceMatrix(m_positions);
    m_parameters     = m_system->getWeights();

    return (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_parameters(m_elementNumber,0) * m_distanceMatrix).cwiseInverse();
}

double PadeJastrow::g(int i, int j, int k, int l) {
    m_positions = m_system->getPositions();
    m_distanceMatrix = calculateDistanceMatrix(m_positions);
    return (m_positions(i) - m_positions(j))/m_distanceMatrix(k,l);
}

double PadeJastrow::beta(int i, int j) {
    m_parameters  = m_system->getWeights();
    int k = m_numberOfDimensions * i + j + 1;      // Stack Gamma matrix
    return m_parameters(m_elementNumber, k);
}

double PadeJastrow::gamma() {
    m_parameters = m_system->getWeights();
    return m_parameters(m_elementNumber, 0);
}

double PadeJastrow::evaluate(Eigen::VectorXd positions) {
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + gamma() * m_distanceMatrix(i,j));
            PadeJastrowFactor += beta(i,j) * f * m_distanceMatrix(i,j);
        }
    }
    return exp(PadeJastrowFactor);
}

double PadeJastrow::evaluateSqrd(Eigen::VectorXd positions) {
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + gamma() * m_distanceMatrix(i,j));
            PadeJastrowFactor += beta(i,j) * f * m_distanceMatrix(i,j);
        }
    }
    return exp(2 * PadeJastrowFactor);
}

double PadeJastrow::computeFirstDerivative(const Eigen::VectorXd positions, int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    Eigen::MatrixXd F = f();

    double derivative = 0;
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += beta(k_p,j_p) * F(k_p,j_p) * F(k_p, j_p) * (positions(k) - positions(j))/m_distanceMatrix(k_p,j_p);
    }

    return derivative;
}

double PadeJastrow::computeSecondDerivative() {
    m_parameters         = m_system->getWeights();
    Eigen::MatrixXd F = f();

    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double G = g(i,j,i_p,j_p);
            derivative += beta(i_p,j_p) * F(i_p, j_p) * F(i_p, j_p) * (1 - G * G * (1 + 3 * gamma() * m_distanceMatrix(i_p,j_p)) * F(i_p, j_p)) / m_distanceMatrix(i_p,j_p);
        }
    }
    return derivative;
}

Eigen::VectorXd PadeJastrow::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    Eigen::MatrixXd F = f();

    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    //Update beta matrix
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        int l = m_numberOfDimensions * k_p + j_p + 1;
        gradients(l) = - 0.5 * F(k_p, j_p) * F(k_p, j_p) * g(k,j,k_p,j_p);
    }

    //Update gamma
    double derivative = 0;
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += beta(k_p,j_p) * F(k_p, j_p) * F(k_p, j_p) * F(k_p, j_p) * (m_positions(k) - m_positions(j));
    }
    gradients(0) = derivative;
    return gradients;
}

Eigen::VectorXd PadeJastrow::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
    Eigen::MatrixXd F = f();

    //Update Beta matrix
    for(int i=0; i < m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double G = g(i,j,i_p,j_p);
            int l = m_numberOfParticles * i_p + j_p + 1;      // Stack Beta matrix
            gradients(l) = 0.5 * F(i_p, j_p) * F(i_p, j_p) * (G * G * (1 + 3 * gamma() * m_distanceMatrix(i_p,j_p)) * F(i_p, j_p) - 1) / m_distanceMatrix(i_p,j_p);
        }
    }

    //Update gamma
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double G = g(i,j,i_p,j_p);
            derivative -= beta(i_p,j_p) * F(i_p, j_p) * F(i_p, j_p) * F(i_p, j_p) * (3 * G * G * F(i_p, j_p) * gamma() * m_distanceMatrix(i_p,j_p) - 1);
        }
    }
    gradients(0) = derivative;
    return gradients;
}
