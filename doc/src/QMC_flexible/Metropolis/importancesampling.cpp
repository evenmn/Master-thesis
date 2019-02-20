#include "importancesampling.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../Math/random2.h"

using std::cout;
using std::endl;

ImportanceSampling::ImportanceSampling(System* system) :
        Metropolis(system) {
    m_numberOfParticles     = m_system->getNumberOfParticles();
    m_numberOfDimensions    = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    //m_diff                  = m_system->getDiffusionConstant();
    m_stepLength            = m_system->getStepLength();
    m_waveFunctionVector    = m_system->getWaveFunctionElements();
}

double ImportanceSampling::QuantumForce(const Eigen::VectorXd positions, int i) {
    double QF = 0;
    for(auto& j : m_waveFunctionVector) {
        QF += j->computeFirstDerivative(positions, i);
    }
    return 2*QF;
}

double ImportanceSampling::GreenFuncSum(const Eigen::VectorXd newPositions) {
    double GreenSum  = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        double GreenFunc = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            double QForceOld = QuantumForce(m_positions, m_numberOfDimensions*i+j);
            double QForceNew = QuantumForce(newPositions, m_numberOfDimensions*i+j);
            GreenFunc += 0.5*(QForceOld + QForceNew) * (0.5*m_diff*m_stepLength*(QForceOld - QForceNew)-newPositions(m_numberOfDimensions*i+j)+m_positions(m_numberOfDimensions*i+j));
        }
        GreenSum += exp(GreenFunc);
    }
    return GreenSum;
}

bool ImportanceSampling::acceptMove() {
    m_positions             = m_system->getParticles();

    Random2 rand;

    int pRand = rand.nextInt(m_numberOfFreeDimensions);

    Eigen::VectorXd newPositions      = m_positions;
    newPositions(pRand) = m_positions(pRand) + m_diff * QuantumForce(m_positions, pRand) * m_stepLength + rand.nextGaussian(0,1) * sqrt(m_stepLength);   //Update position                                 //Update v

    double psiOld = m_system->evaluateWaveFunctionSqrd();
    m_system->updateAllArrays(newPositions, pRand);
    double psiNew = m_system->evaluateWaveFunctionSqrd();

    double w = GreenFuncSum(newPositions) * (psiNew/psiOld);
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
