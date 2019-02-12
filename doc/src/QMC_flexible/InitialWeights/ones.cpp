#include "ones.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

Ones::Ones(System* system)  :  InitialWeights(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfElements = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    setupInitialWeights();
}

void Ones::setupInitialWeights() {
    m_parameters = 2*Eigen::MatrixXd::Ones(m_numberOfElements, m_maxNumberOfParametersPerElement);
}
