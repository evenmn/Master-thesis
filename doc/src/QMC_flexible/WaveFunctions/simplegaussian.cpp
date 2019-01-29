#include "simplegaussian.h"
#include <cassert>
#include <iostream>
#include "wavefunction.h"
#include "../system.h"

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(Eigen::MatrixXd particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrdElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrdElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrdElementWise;
    }

    return exp(-0.5 * m_parameters.at(0) * r.sum());
}

double SimpleGaussian::computeFirstDerivative(Eigen::MatrixXd particles, int k) {
    // Calculating the kinetic energy term, -0.5 * laplacian
    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }

    return -m_parameters.at(0) * r(k);
}

double SimpleGaussian::computeSecondDerivative(Eigen::MatrixXd particles) {
    return -m_parameters.at(0) * particles.rows() * particles.cols();
}

double SimpleGaussian::computeFirstEnergyDerivative(Eigen::MatrixXd particles) {
    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }
    return 0.5 * m_numberOfParticles * m_numberOfDimensions -m_parameters.at(0) * r.sum() * r.sum();
}

double SimpleGaussian::computeSecondEnergyDerivative(Eigen::MatrixXd particles) {
    return 0.5 * particles.rows() * particles.cols();
}
