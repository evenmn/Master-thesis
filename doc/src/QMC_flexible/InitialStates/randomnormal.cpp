#include "randomnormal.h"
#include <iostream>
#include <cassert>
#include "Math/random.h"
#include "../system.h"
#include "WaveFunctions/wavefunction.h"

RandomNormal::RandomNormal(System*    system)  :
        InitialState(system) {
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

void RandomNormal::setupInitialState() {
    Random rand;
    m_positions = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for (int i=0; i < m_numberOfFreeDimensions; i++) {
        m_positions(i) = rand.nextGaussian(0,1);
    }
    for(auto& i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions);
    }
}
