#include "atomicnucleus.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

AtomicNucleus::AtomicNucleus(System* system, int Z) :
        Hamiltonian(system) {
    assert(Z > 0);
    m_Z                     = Z;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
    m_particles             = m_system->getParticles();
}

double AtomicNucleus::computeLocalEnergy() {
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
            sqrtElementWise += m_particles(i,j) * m_particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
        r_inv(i) = 1/r(i);
    }

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double sqrtElementWise = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                double numb = m_particles(i,d) - m_particles(j,d);
                sqrtElementWise += numb * numb;
            }
            R(i,j) = sqrt(sqrtElementWise);
        }
    }

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=0; j<i; j++) {
                interactionEnergy += 1/R(i,j);
            }
        }
    }

    double externalEnergy = -0.5 * m_Z * r_inv.sum();
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
