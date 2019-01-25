#include "particle.h"
#include <cassert>

Particle::Particle() {
}

//void Particle::setPosition(const std::vector<double> &position) {
void Particle::setPosition(const Eigen::MatrixXd &position) {
    assert(position.size() == m_numberOfDimensions);
    m_position = position;
}

void Particle::adjustPosition(double change, int dimension) {
    //m_position.at(dimension) += change;
    m_position(dimension,dimension) += change;
}

void Particle::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}
