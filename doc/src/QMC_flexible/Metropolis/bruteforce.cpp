#include "bruteforce.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../Math/random2.h"

using std::cout;
using std::endl;

BruteForce::BruteForce(System* system) :
        Metropolis(system) {
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_stepLength            = m_system->getStepLength();
}

bool BruteForce::acceptMove() {
    m_positions             = m_system->getParticles();


    Random2 rand;

    int pRand = rand.nextInt(m_numberOfFreeDimensions);

    Eigen::VectorXd newPositions      = m_positions;

    newPositions(pRand) = m_positions(pRand) + (rand.nextDouble() - 0.5) * m_stepLength;

    double psiOld = m_system->evaluateWaveFunctionSqrd();
    m_system->updateAllArrays(newPositions, pRand);
    double psiNew = m_system->evaluateWaveFunctionSqrd();

    double w = psiNew/psiOld;
    double r = rand.nextDouble();
    if(w > r) {
        m_positions(pRand)        = newPositions(pRand);
        return true;
    }
    else {
        m_system->resetAllArrays();
        return false;
    }
}
