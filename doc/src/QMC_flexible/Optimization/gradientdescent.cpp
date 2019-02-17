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
    m_numberOfFreeDimensions          = m_system->getNumberOfFreeDimensions();
    m_numberOfWaveFunctionElements    = m_system->getNumberOfWaveFunctionElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxNumberOfParametersPerElement();
    m_waveFunctionVector              = m_system->getWaveFunctionElements();
    m_eta                             = m_system->getLearningRate();
}

Eigen::VectorXd GradientDescent::getImmediateGradients(WaveFunction* waveFunction) {
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

Eigen::MatrixXd GradientDescent::getAllImmediateGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    for(int i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(i) += getImmediateGradients(m_waveFunctionVector[unsigned(i)]);
    }


    Eigen::MatrixXd gradients2 = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    m_positions = m_system->getPositions();
    Eigen::MatrixXd parameters = m_system->getWeights();

    gradients2.row(0).head(4) = m_positions.transpose() - parameters.row(0).head(4);

    Eigen::VectorXd e = Eigen::VectorXd::Zero(2);
    for(int i=0; i<2; i++) {
        double Sum = 0;
        for(int j=0; j<4; j++) {
            Sum += parameters(1, i*2 + j) * m_positions(j);
        }
        e(i) = 1/(1 + Sum);
    }

    gradients2.row(1).head(2) = e;

    for(int i=0; i<8; i++) {
        gradients2(1, i+2) = (m_positions * e.transpose())(i%4, int(i/4));
    }

    //cout << gradients << endl;
    //cout << endl;
    //cout << gradients2 << endl;
    //cout << endl;
    //cout << endl;

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
