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

    double interactionEnergy = 0;
    if(m_interaction) {
        for(int i=0; i<m_numberOfParticles; i++) {
            for(int j=0; j<i; j++) {
                double sqrdElementWise = 0;
                for(int d=0; d<m_numberOfDimensions; d++) {
                    double numb = m_positions(i*m_numberOfDimensions + d) - m_positions(j*m_numberOfDimensions + d);
                    sqrdElementWise += numb * numb;
                }
                interactionEnergy += 1/sqrt(sqrdElementWise);
            }
        }
    }

    double nucleusEnergy = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int d=0; d<m_numberOfDimensions; d++) {
            sqrtElementWise += m_positions(i*m_numberOfDimensions + d) * m_positions(i*m_numberOfDimensions + d);
        }
        nucleusEnergy += 1/sqrt(sqrtElementWise);
    }


    /*
    double interactionEnergy = 0;
    if(m_interaction) {
        Eigen::MatrixXd Inverse = m_distanceMatrix.cwiseInverse().unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; });
        interactionEnergy = Inverse.sum();
    }
    */

    double externalEnergy = - 0.5 * m_Z * nucleusEnergy;
    double kineticEnergy  = m_system->getKineticEnergy();

    return kineticEnergy + externalEnergy + interactionEnergy;
}
