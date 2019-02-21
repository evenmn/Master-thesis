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

double PadeJastrow::calculateDistanceMatrixElement(int i, int j, Eigen::VectorXd particles) {
    // Update element (i,j) in Dist_inv_invance matrix

    double dist = 0;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = particles(m_numberOfDimensions*i+d)-particles(m_numberOfDimensions*j+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

Eigen::MatrixXd PadeJastrow::calculateDistanceMatrix(Eigen::VectorXd particles) {
    bool loops = true;

    if(loops) {
        Eigen::MatrixXd distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=0; j<i; j++) {
                distanceMatrix(i,j) = calculateDistanceMatrixElement(i,j,particles);
            }
        }
        return distanceMatrix;
    }
    else if((loops=false)) {
        Eigen::Map<Eigen::MatrixXd> r(particles.data(), m_numberOfParticles, m_numberOfDimensions);
        return Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    }
}

void PadeJastrow::calculateDistanceMatrixCross(int par, Eigen::VectorXd particles, Eigen::MatrixXd &distanceMatrix) {
    // Update distance matrix when position of particle "par" is changed

    // Update row
    for(int i=0; i<par; i++) {
        distanceMatrix(par, i) = calculateDistanceMatrixElement(par, i, particles);
    }
    // Update column
    for(int i=par+1; i<m_numberOfParticles; i++) {
        distanceMatrix(i, par) = calculateDistanceMatrixElement(i, par, particles);
    }
}

void PadeJastrow::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_distanceMatrix = calculateDistanceMatrix(positions);
    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();
    m_g     = Eigen::MatrixXd::Zero(m_numberOfFreeDimensions, m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = (m_positions(i) - m_positions(j))/m_distanceMatrix(int(i/m_numberOfDimensions),int(j/m_numberOfDimensions));
        }
    }
}

void PadeJastrow::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;

    m_oldDistanceMatrix = m_distanceMatrix;
    calculateDistanceMatrixCross(int(pRand/m_numberOfDimensions), positions, m_distanceMatrix);

    m_oldF  = m_f;
    m_f     = (Eigen::MatrixXd::Ones(m_numberOfParticles, m_numberOfParticles) + m_gamma * m_distanceMatrix).cwiseInverse();

    m_oldG  = m_g;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        for(int j=0; j<m_numberOfFreeDimensions; j++) {
            m_g(i,j) = (m_positions(i) - m_positions(j))/m_distanceMatrix(int(i/m_numberOfDimensions),int(j/m_numberOfDimensions));
        }
    }
}

void PadeJastrow::resetArrays() {
    m_positions = m_oldPositions;
    m_distanceMatrix = m_oldDistanceMatrix;
}

void PadeJastrow::updateParameters(Eigen::MatrixXd parameters) {
    Eigen::VectorXd X = parameters.row(m_elementNumber).segment(1, m_numberOfParticles*m_numberOfParticles);
    Eigen::Map<Eigen::MatrixXd> beta(X.data(), m_numberOfParticles, m_numberOfParticles);
    m_beta  = beta.transpose();
    m_gamma = parameters(m_elementNumber, 0);
}

double PadeJastrow::evaluate() {
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_gamma * m_distanceMatrix(i,j));
            PadeJastrowFactor += m_beta(i,j) * f * m_distanceMatrix(i,j);
        }
    }
    return exp(PadeJastrowFactor);
}

double PadeJastrow::evaluateSqrd() {
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_gamma * m_distanceMatrix(i,j));
            PadeJastrowFactor += m_beta(i,j) * f * m_distanceMatrix(i,j);
        }
    }
    return exp(2 * PadeJastrowFactor);
}

double PadeJastrow::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    double derivative = 0;
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p,j_p) * m_f(k_p,j_p) * m_f(k_p, j_p) * (positions(k) - positions(j))/m_distanceMatrix(k_p,j_p);
    }

    return derivative;
}

double PadeJastrow::computeSecondDerivative() {
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double G = m_g(i,j);
            derivative += m_beta(i_p,j_p) * m_f(i_p, j_p) * m_f(i_p, j_p) * (1 - G * G * (1 + 3 * m_gamma * m_distanceMatrix(i_p,j_p)) * m_f(i_p, j_p)) / m_distanceMatrix(i_p,j_p);
        }
    }
    return derivative;
}

Eigen::VectorXd PadeJastrow::computeFirstEnergyDerivative(int k) {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    int k_p = int(k/m_numberOfDimensions);  //Particle associated with k
    int k_d = k%m_numberOfDimensions;       //Dimension associated with k

    //Update beta matrix
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        int l = m_numberOfDimensions * k_p + j_p + 1;
        gradients(l) = - 0.5 * m_f(k_p, j_p) * m_f(k_p, j_p) * m_g(k,j);
    }

    //Update gamma
    double derivative = 0;
    for(int j_p=0; j_p<k_p; j_p++) {
        int j = j_p * m_numberOfDimensions + k_d;
        derivative += m_beta(k_p,j_p) * m_f(k_p, j_p) * m_f(k_p, j_p) * m_f(k_p, j_p) * (m_positions(k) - m_positions(j));
    }
    gradients(0) = derivative;
    return gradients;
}

Eigen::VectorXd PadeJastrow::computeSecondEnergyDerivative() {
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    //Update Beta matrix
    for(int i=0; i < m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double G = m_g(i,j);
            int l = m_numberOfParticles * i_p + j_p + 1;      // Stack Beta matrix
            gradients(l) = 0.5 * m_f(i_p, j_p) * m_f(i_p, j_p) * (G * G * (1 + 3 * m_gamma * m_distanceMatrix(i_p,j_p)) * m_f(i_p, j_p) - 1) / m_distanceMatrix(i_p,j_p);
        }
    }

    //Update gamma
    double derivative = 0;
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        int i_p = int(i/m_numberOfDimensions);  //Particle associated with k
        int i_d = i%m_numberOfDimensions;       //Dimension associated with k
        for(int j_p=0; j_p<i_p; j_p++) {
            int j = j_p * m_numberOfDimensions + i_d;
            double G = m_g(i,j);
            derivative -= m_beta(i_p,j_p) * m_f(i_p, j_p) * m_f(i_p, j_p) * m_f(i_p, j_p) * (3 * G * G * m_f(i_p, j_p) * m_gamma * m_distanceMatrix(i_p,j_p) - 1);
        }
    }
    gradients(0) = derivative;
    return gradients;
}
