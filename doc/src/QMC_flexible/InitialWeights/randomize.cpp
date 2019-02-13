#include "randomize.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

Randomize::Randomize(System*    system)  :  InitialWeights(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfElements   = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    setupInitialWeights();
}

void Randomize::setupInitialWeights() {
    m_parameters = Eigen::MatrixXd::Random(m_numberOfElements, m_maxNumberOfParametersPerElement);
}
