#include "randomize.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

Randomize::Randomize(System*    system, int numberOfElements)  :  InitialWeights(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfElements   = numberOfElements;
    setupInitialWeights();
}

void Randomize::setupInitialWeights() {
    int maxNumberOfParametersPerElement = m_numberOfParticles * m_numberOfParticles + m_numberOfParticles;
    m_parameters = Eigen::MatrixXd::Random(m_numberOfElements, maxNumberOfParametersPerElement);
}
