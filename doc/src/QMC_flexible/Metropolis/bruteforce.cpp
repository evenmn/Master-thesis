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
    m_diff                  = m_system->getDiffusionConstant();
    m_stepLength            = m_system->getStepLength();
}

bool BruteForce::acceptMove() {
    m_particles             = m_system->getParticles();
    m_radialVector          = m_system->getRadialVector();
    m_distanceMatrix        = m_system->getDistanceMatrix();

    Random2 rand;

    int pRand = rand.nextInt(m_numberOfFreeDimensions);

    Eigen::VectorXd newPositions      = m_particles;
    Eigen::VectorXd newRadialVector   = m_radialVector;
    Eigen::MatrixXd newDistanceMatrix = m_distanceMatrix;

    newPositions(pRand) = m_particles(pRand) + (rand.nextDouble() - 0.5) * m_stepLength;
    //calculateDistanceMatrixCross(int(pRand/m_numberOfDimensions), newPositions, newDistanceMatrix);

    newRadialVector   = m_system->calculateRadialVector(newPositions);
    newDistanceMatrix = m_system->calculateDistanceMatrix(newPositions);

    double psiOld = m_system->evaluateWaveFunctionSqrd(m_particles, m_radialVector, m_distanceMatrix);
    double psiNew = m_system->evaluateWaveFunctionSqrd(newPositions, newRadialVector, newDistanceMatrix);

    double w = psiNew/psiOld;
    double r = rand.nextDouble();

    if(w > r) {
        m_particles(pRand)        = newPositions(pRand);
        m_radialVector            = newRadialVector;
        m_distanceMatrix          = newDistanceMatrix;
        return true;
    }
    return false;
}
