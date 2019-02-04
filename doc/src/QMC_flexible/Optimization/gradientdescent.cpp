#include "gradientdescent.h"
#include <cassert>
#include "../system.h"
//#include "../WaveFunctions/wavefunction.h"

GradientDescent::GradientDescent(System* system) :
        Optimization(system) {
}

void GradientDescent::getEnergyGradient(double EL_avg, Eigen::MatrixXd grad_tot, Eigen::MatrixXd gradE_tot, Eigen::MatrixXd &gradients) {
    gradients = 2 * (gradE_tot - EL_avg * grad_tot)/m_system->getNumberOfMetropolisSteps();
}
