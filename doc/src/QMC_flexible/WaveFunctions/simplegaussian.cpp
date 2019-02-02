#include "simplegaussian.h"
#include <cassert>
#include <iostream>
//#include "wavefunction.h"
#include "../system.h"

SimpleGaussian::SimpleGaussian(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber      = elementNumber;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

double SimpleGaussian::evaluate(Eigen::MatrixXd particles) {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrdElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrdElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrdElementWise);
    }
    return exp(-0.5 * m_parameters(m_elementNumber, 0) * (r.cwiseAbs2()).sum());
}

double SimpleGaussian::computeFirstDerivative(int k) {
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

    return -m_parameters(m_elementNumber, 0) * r(k);
}

double SimpleGaussian::computeSecondDerivative() {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    return -m_parameters(m_elementNumber, 0) * m_numberOfParticles * m_numberOfDimensions;
}

double SimpleGaussian::computeFirstEnergyDerivative() {
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
    return 0.5 * m_numberOfParticles * m_numberOfDimensions - m_parameters(m_elementNumber, 0) * r.sum() * r.sum();
}

double SimpleGaussian::computeSecondEnergyDerivative() {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    return 0.5 * m_numberOfParticles * m_numberOfDimensions;
}
