#include "randomnormal.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

RandomNormal::RandomNormal(System*    system)  :
        InitialState(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    setupInitialState();
}

void RandomNormal::setupInitialState() {
    Random rand;
    Eigen::MatrixXd positions = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfDimensions);
    for (int i=0; i < m_numberOfParticles; i++) {
        for (int j=0; j < m_numberOfDimensions; j++) {
            positions(i,j) = rand.nextGaussian(0,1);
        }
    }
    m_particles = positions;

    m_radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    m_system->calculateRadialVector(m_particles, m_radialVector);
    m_system->calculateDistanceMatrix(m_particles, m_distanceMatrix);
}
