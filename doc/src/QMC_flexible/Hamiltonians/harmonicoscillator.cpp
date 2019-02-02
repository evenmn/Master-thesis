#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega                 = omega;
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
    m_particles             = m_system->getParticles();
}

double HarmonicOscillator::computeLocalEnergy() {
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    std::cout << m_particles << std::endl;

    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            std::cout << "gkgkgk" << std::endl;
            sqrtElementWise += m_particles(i,j) * m_particles(i,j);
            std::cout << "gkgkgk" << std::endl;
        }
        r(i) = sqrt(sqrtElementWise);
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

    double externalEnergy = 0.5 * m_omega * m_omega * (r.cwiseAbs2()).sum();
    double kineticEnergy  = m_system->getKineticEnergy();

    /*
    std::cout << kineticEnergy << std::endl;
    std::cout << externalEnergy << std::endl;
    std::cout << interactionEnergy << std::endl;
    std::cout << " " << std::endl;
    */

    return kineticEnergy + externalEnergy + interactionEnergy;
}

