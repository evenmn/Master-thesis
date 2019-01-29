#include "gradientdescent.h"
#include <cassert>
#include "optimization.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"

GradientDescent::GradientDescent(System* system) :
        Optimization(system) {
    //assert(alpha >= 0);
    //m_numberOfParameters = 1;
    //m_parameters.reserve(1);
    //m_parameters.push_back(alpha);
}

double GradientDescent::gradient(Eigen::MatrixXd particles) {
    double gradient = 0; //m_system->getWaveFunction()->computeFirstEnergyDerivative(particles);

    return gradient;
}
