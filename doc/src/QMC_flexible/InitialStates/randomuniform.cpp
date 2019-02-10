#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

RandomUniform::RandomUniform(System*    system)  :
        InitialState(system) {
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    Random rand;
    Eigen::VectorXd positions = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for (int i=0; i < m_numberOfFreeDimensions; i++) {
        positions(i) = rand.nextDouble();
    }
    m_particles = positions;

    m_distanceMatrix = m_system->calculateDistanceMatrix(m_particles);
    m_radialVector   = m_system->calculateRadialVector(m_particles);
}
