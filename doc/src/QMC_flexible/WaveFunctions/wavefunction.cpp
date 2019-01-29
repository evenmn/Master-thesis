#include "wavefunction.h"
#include <iostream>

WaveFunction::WaveFunction(std::vector<System *> system) {
    m_system = system;
}

double WaveFunction::TotalEvaluation(Eigen::MatrixXd particles) {
    double TotalWF = 0;
    for (unsigned i=0; i < m_system.size(); i++) {
        std::cout << m_system[i] << std::endl;
        //TotalWF += m_system[i]->evaluate(particles);
    }
    return TotalWF;
}

WaveFunction::~WaveFunction() {};
