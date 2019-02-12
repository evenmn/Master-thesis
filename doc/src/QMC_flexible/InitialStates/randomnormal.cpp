#include "randomnormal.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

RandomNormal::RandomNormal(System*    system)  :
        InitialState(system) {
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

void RandomNormal::setupInitialState() {
    Random rand;
    Eigen::VectorXd positions = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for (int i=0; i < m_numberOfFreeDimensions; i++) {
        positions(i) = rand.nextGaussian(0,1);
    }
    m_positions = positions;
    m_distanceMatrix = m_system->calculateDistanceMatrix(m_positions);
    m_radialVector   = m_system->calculateRadialVector(m_positions);
}
