#include "ones.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"

Ones::Ones(System*    system,
                             unsigned numberOfElements)  :
        InitialWeights(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfElements   = numberOfElements;
    setupInitialWeights();
}

void Ones::setupInitialWeights() {
    int maxNumberOfParametersPerElement = m_numberOfParticles*m_numberOfParticles + m_numberOfParticles;
    m_parameters = Eigen::MatrixXd::Ones(m_numberOfElements, maxNumberOfParametersPerElement);
}
