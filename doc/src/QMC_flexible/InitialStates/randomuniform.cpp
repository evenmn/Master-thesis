#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

RandomUniform::RandomUniform(System*    system)  :
        InitialState(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    Random rand;
    Eigen::MatrixXd positions = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfDimensions);
    for(int i=0; i < m_numberOfParticles; i++) {
        for(int j=0; j < m_numberOfDimensions; j++) {
            positions(i,j) = rand.nextDouble();
        }
    }
    m_particles = positions;

    m_radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    m_system->calculateRadialVector(m_particles, m_radialVector);
    m_system->calculateDistanceMatrix(m_particles, m_distanceMatrix);
}
