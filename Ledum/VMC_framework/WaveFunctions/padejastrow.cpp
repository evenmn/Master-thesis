#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

PadeJastrow::PadeJastrow(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double PadeJastrow::evaluate(Eigen::MatrixXd particles) {
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
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }

    return exp(-m_parameters.at(0) * (r.cwiseAbs2()).sum());
}

double PadeJastrow::computeFirstDerivative(Eigen::MatrixXd particles, int k) {
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

double PadeJastrow::computeDoubleDerivative(Eigen::MatrixXd particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    return -m_parameters.at(0) * particles.rows() * particles.cols();
}

double PadeJastrow::computeFirstEnergyDerivative(Eigen::MatrixXd particles) {
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
    return 0.5 * double(r.transpose()*r);
}

double PadeJastrow::computeDoubleEnergyDerivative(Eigen::MatrixXd particles) {
    return 0.5 * particles.rows() * particles.cols();
}
