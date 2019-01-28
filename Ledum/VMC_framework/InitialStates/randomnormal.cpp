#include "randomnormal.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

RandomNormal::RandomNormal(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomNormal::setupInitialState() {
    Random rand;
    Eigen::MatrixXd positions = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfDimensions);
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfDimensions; j++) {
            positions(i,j) = rand.nextGaussian(0,1);
        }
        m_particles = positions;
    }
    //std::cout << m_particles << std::endl;
}
