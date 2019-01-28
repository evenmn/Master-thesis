#include "hydrogenorbital.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

HydrogenOrbital::HydrogenOrbital(System* system, double beta) :
        WaveFunction(system) {
    assert(beta >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(beta);
}

double HydrogenOrbital::evaluate(Eigen::MatrixXd particles) {
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

    return exp(-m_parameters.at(0)* m_numberOfParticles * r.sum());
}

double HydrogenOrbital::computeFirstDerivative(Eigen::MatrixXd particles, int k) {

    long m_numberOfParticles = particles.rows();
    return -m_parameters.at(0) * m_numberOfParticles;
}

double HydrogenOrbital::computeDoubleDerivative(Eigen::MatrixXd particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */
    return 0;
}

double HydrogenOrbital::computeFirstEnergyDerivative(Eigen::MatrixXd particles) {
    long m_numberOfParticles = particles.rows();
    return 0.5 * m_numberOfParticles * m_numberOfParticles;
}

double HydrogenOrbital::computeDoubleEnergyDerivative(Eigen::MatrixXd particles) {
    return 0;
}
