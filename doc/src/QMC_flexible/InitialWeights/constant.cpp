#include "constant.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

Constant::Constant(System* system, double factor)  :  InitialWeights(system) {
    m_numberOfDimensions              = m_system->getNumberOfDimensions();
    m_numberOfParticles               = m_system->getNumberOfParticles();
    m_numberOfElements                = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_factor                          = factor;
    setupInitialWeights();
}

void Constant::setupInitialWeights() {
    m_parameters = m_factor * Eigen::MatrixXd::Ones(m_numberOfElements, m_maxNumberOfParametersPerElement);
    m_system->updateAllParameters(m_parameters);
}
