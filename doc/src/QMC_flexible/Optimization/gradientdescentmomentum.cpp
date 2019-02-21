#include "gradientdescentmomentum.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../Math/random2.h"

GradientDescentMomentum::GradientDescentMomentum(System* system, double gamma) :
        Optimization(system) {
    m_numberOfFreeDimensions          = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
    m_v                               = Eigen::MatrixXd::Ones(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_gamma                           = gamma;
}

Eigen::VectorXd GradientDescentMomentum::getImmediateGradients(WaveFunction* waveFunction) {
    m_positions = m_system->getParticles();
    Eigen::VectorXd TotalGradients = waveFunction->computeSecondEnergyDerivative();
    for(int k = 0; k < m_numberOfFreeDimensions; k++) {
        double Sum = 0;
        for(auto& i : m_waveFunctionVector) {
            Sum += i->computeFirstDerivative(m_positions, k);
        }
        TotalGradients += 2 * Sum * waveFunction->computeFirstEnergyDerivative(k);
    }
    return TotalGradients;
}

Eigen::MatrixXd GradientDescentMomentum::getAllImmediateGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    for(int i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(i) += getImmediateGradients(m_waveFunctionVector[unsigned(i)]);
    }
    return gradients;
}

Eigen::MatrixXd GradientDescentMomentum::getEnergyGradient() {
    double E            = m_system->getSampler()->getEnergy();
    Eigen::MatrixXd dE  = m_system->getSampler()->getdE();
    Eigen::MatrixXd dEE = m_system->getSampler()->getdEE();
    return 2 * (dEE - E * dE)/ int((1 - m_system->getEquilibrationFraction()) * m_system->getNumberOfMetropolisSteps());
}

Eigen::MatrixXd GradientDescentMomentum::updateParameters() {
    m_v = m_gamma * m_v + m_eta * getEnergyGradient();
    return m_v;
}
