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
}

double AtomicNucleus::computeLocalEnergy() {
    m_positions             = m_system->getParticles();
    m_radialVector          = m_system->getRadialVector();
    m_distanceMatrix        = m_system->getDistanceMatrix();

    double radialVectorInverseSum = 0;
    for(int i=0; i < m_numberOfParticles; i++) {
        radialVectorInverseSum += 1/m_radialVector(i);
    }

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=0; j<i; j++) {
                interactionEnergy += 1/m_distanceMatrix(i,j);
            }
        }
    }

    double externalEnergy = -0.5 * m_Z * radialVectorInverseSum;
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
