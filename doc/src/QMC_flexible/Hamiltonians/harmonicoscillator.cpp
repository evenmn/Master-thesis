#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system) :
        Hamiltonian(system) {
    m_omega                 = m_system->getFrequency();
    assert(m_omega > 0);
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_interaction           = m_system->getInteraction();
}

double HarmonicOscillator::computeLocalEnergy() {
    m_particles             = m_system->getParticles();
    m_radialVector          = m_system->getRadialVector();
    m_distanceMatrix        = m_system->getDistanceMatrix();

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=0; j<i; j++) {
                interactionEnergy += 1/m_distanceMatrix(i,j);
            }
        }
    }

    double externalEnergy = 0.5 * m_omega * m_omega * (m_particles.cwiseAbs2()).sum();
    double kineticEnergy  = m_system->getKineticEnergy();

    /*
    std::cout << kineticEnergy << std::endl;
    std::cout << externalEnergy << std::endl;
    std::cout << interactionEnergy << std::endl;
    std::cout << " " << std::endl;
    */

    return kineticEnergy + externalEnergy + interactionEnergy;
}

