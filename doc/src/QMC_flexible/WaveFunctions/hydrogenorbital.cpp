#include "hydrogenorbital.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

HydrogenOrbital::HydrogenOrbital(System* system, double beta) :
        WaveFunction(system) {
    assert(beta >= 0);
    m_numberOfParameters = 1;
    //m_parameters.reserve(1);
    //m_parameters.push_back(beta);
    m_beta = beta;
}

Eigen::VectorXd HydrogenOrbital::calculateRadialVector(Eigen::VectorXd particles) {
    Eigen::VectorXd radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int d=0; d<m_numberOfDimensions; d++) {
            sqrtElementWise += particles(i*m_numberOfDimensions + d) * particles(i*m_numberOfDimensions + d);
        }
        radialVector(i) = sqrt(sqrtElementWise);
    }
    return radialVector;
}

void HydrogenOrbital::initializeArrays(Eigen::VectorXd positions) {
    m_positions       = positions;
    m_radialVector    = calculateRadialVector(positions);
}

void HydrogenOrbital::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions    = m_positions;
    m_positions       = positions;

    m_oldRadialVector = m_radialVector;
    m_radialVector    = calculateRadialVector(positions);
}

void HydrogenOrbital::resetArrays() {
    m_positions       = m_oldPositions;
    m_radialVector    = m_oldRadialVector;
}

void HydrogenOrbital::updateParameters(Eigen::MatrixXd parameters) {
    //m_a = (parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
}

double HydrogenOrbital::evaluate(Eigen::MatrixXd positions) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    long m_numberOfParticles = positions.rows();
    long m_numberOfDimensions = positions.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += positions(i,j) * positions(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }

    return exp(-m_beta* m_numberOfParticles * r.sum());
}

double HydrogenOrbital::computeDerivative(Eigen::MatrixXd positions) {
    // Calculating the kinetic energy term, -0.5 * laplacian
    return -0.5 * m_beta * m_beta;
}

double HydrogenOrbital::computeEnergyDerivative(Eigen::MatrixXd positions) {
    return -m_beta;
}
