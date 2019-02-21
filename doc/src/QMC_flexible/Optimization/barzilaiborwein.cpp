#include "barzilaiborwein.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../sampler.h"
#include "../WaveFunctions/wavefunction.h"
#include "../InitialWeights/initialweights.h"
#include "../Math/random2.h"

BarzilaiBorwein::BarzilaiBorwein(System* system) :
        Optimization(system) {
    m_numberOfFreeDimensions          = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
    m_parameters                      = m_system->getInitialWeights()->getWeights();
    m_gradients                       = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_oldParameters                   = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
}

Eigen::VectorXd BarzilaiBorwein::getImmediateGradients(WaveFunction* waveFunction) {
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

Eigen::MatrixXd BarzilaiBorwein::getAllImmediateGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    for(int i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(i) += getImmediateGradients(m_waveFunctionVector[unsigned(i)]);
    }
    return gradients;
}

Eigen::MatrixXd BarzilaiBorwein::getEnergyGradient() {
    double E            = m_system->getSampler()->getEnergy();
    Eigen::MatrixXd dE  = m_system->getSampler()->getdE();
    Eigen::MatrixXd dEE = m_system->getSampler()->getdEE();
    return 2 * (dEE - E * dE)/ int((1 - m_system->getEquilibrationFraction()) * m_system->getNumberOfMetropolisSteps());
}

Eigen::MatrixXd BarzilaiBorwein::updateParameters() {
    m_oldGradients = m_gradients;
    m_parameters    = m_system->getWeights();
    m_gradients     = getEnergyGradient();
    Eigen::MatrixXd learningRate = (m_parameters - m_oldParameters).cwiseProduct((m_gradients - m_oldGradients).cwiseInverse());
    for(int i=0; i<m_numberOfWaveFunctionElements; i++) {
        for(int j=0; j<m_maxNumberOfParametersPerElement; j++) {
            if(std::isnan(learningRate(i,j)) || std::isinf(learningRate(i,j))) {
                learningRate(i,j) = m_eta;
            }
        }
    }
    m_oldParameters = m_parameters;
    return m_eta * learningRate.cwiseProduct(m_gradients);
}
