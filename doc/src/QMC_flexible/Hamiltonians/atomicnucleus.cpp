#include "atomicnucleus.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

AtomicNucleus::AtomicNucleus(System* system, double omega, int numberOfParticles, int numberOfDimensions) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
}

double AtomicNucleus::computeLocalEnergy(Eigen::MatrixXd particles) {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    //std::cout << particles << std::endl;

    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    Eigen::VectorXd r_inv = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
        r_inv(i) = 1/r(i);
    }

    double interactionEnergy = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            interactionEnergy += 1/fabs(r(i) - r(j));
        }
    }
    double externalEnergy = -0.5 * m_numberOfParticles * r_inv.sum();
    double kineticEnergy  = 0; //m_system->getWaveFunction()->computeSecondDerivative(particles);
    return -0.5 * kineticEnergy + externalEnergy + interactionEnergy;
}
