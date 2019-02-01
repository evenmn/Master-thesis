#include "ones.h"
#include <cassert>
#include "../system.h"

Ones::Ones(System*    system, unsigned numberOfElements)  :
        InitialWeights(system) {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfElements   = numberOfElements;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    setupInitialState();
}

void Ones::setupInitialState() {
    Eigen::MatrixXd parameters = Eigen::MatrixXd::Ones(m_numberOfElements, 2*m_numberOfParticles*m_numberOfParticles);
    m_parameters = parameters;
}
