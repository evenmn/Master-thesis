#include "gradientdescent.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../Math/random2.h"

using std::cout;
using std::endl;

GradientDescent::GradientDescent(System* system) :
        Optimization(system) {
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector = m_system->getWaveFunctionElements();
    m_eta                    = m_system->getLearningRate();
}

Eigen::VectorXd GradientDescent::getImmediateGradients(WaveFunction* waveFunction) {
    m_positions = m_system->getParticles();
    Eigen::VectorXd TotalGradients = waveFunction->computeSecondEnergyDerivative();
    for(int k = 0; k < m_numberOfFreeDimensions; k++) {
        TotalGradients += 2 * waveFunction->computeFirstDerivative(m_positions, k) * waveFunction->computeFirstEnergyDerivative(k);
    }
    return TotalGradients;
}

Eigen::MatrixXd GradientDescent::getAllImmediateGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    for(int i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(i) += getImmediateGradients(m_waveFunctionVector[unsigned(i)]);
    }
    return gradients;
}

Eigen::MatrixXd GradientDescent::getEnergyGradient() {
    double E            = m_system->getSampler()->getEnergy();
    Eigen::MatrixXd dE  = m_system->getSampler()->getdE();
    Eigen::MatrixXd dEE = m_system->getSampler()->getdEE();
    return 2 * (dEE - E * dE)/ int((1 - m_system->getEquilibrationFraction()) * m_system->getNumberOfMetropolisSteps());
}

Eigen::MatrixXd GradientDescent::updateParameters() {
    return m_eta * getEnergyGradient();
}
