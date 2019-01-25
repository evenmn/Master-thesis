#include "particle.h"
#include <cassert>

Particle::Particle() {
}

void Particle::setPosition(const Eigen::VectorXd &position) {
    assert(position.size() == m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    m_position(dimension) += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
