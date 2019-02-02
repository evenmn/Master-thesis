#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

PadeJastrow::PadeJastrow(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber      = elementNumber;
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
}

double PadeJastrow::evaluate(Eigen::MatrixXd particles) {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);

    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double sqrtElementWise = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                double numb = particles(i,d) - particles(j,d);
                sqrtElementWise += numb * numb;
            }
            R(i,j) = sqrt(sqrtElementWise);
        }
    }

    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            int l = m_numberOfParticles*i + j + 1;      // Stack Gamma matrix
            PadeJastrowFactor += m_parameters(m_elementNumber, l) * R(i,j)/(1 + m_parameters(m_elementNumber, 0) * R(i,j));
        }
    }

    return exp(PadeJastrowFactor);
}

double PadeJastrow::computeFirstDerivative(int k) {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += m_particles(i,j) * m_particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double sqrtElementWise = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                double numb = m_particles(i,d) - m_particles(j,d);
                sqrtElementWise += numb * numb;
            }
            R(i,j) = sqrt(sqrtElementWise);
        }
    }

    double derivative = 0;
    for(int j=0; j<k; j++) {
        double f = 1/(1 + m_parameters(m_elementNumber, 0) * R(k,j));
        int l = m_numberOfParticles*k + j + 1;      // Stack Gamma matrix
        derivative+= m_parameters(m_elementNumber, l) * f * f;
    }
    return derivative;
}

double PadeJastrow::computeSecondDerivative() {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += m_particles(i,j) * m_particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double sqrtElementWise = 0;
            for(int d=0; d<m_numberOfDimensions; d++) {
                double numb = m_particles(i,d) - m_particles(j,d);
                sqrtElementWise += numb * numb;
            }
            R(i,j) = sqrt(sqrtElementWise);
        }
    }

    double derivative = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_parameters(m_elementNumber, 0) * R(i,j));
            int l = m_numberOfParticles*i + j + 1;      // Stack Gamma matrix
            derivative += m_parameters(m_elementNumber, l) * f * f * (1/r(i) - m_parameters(m_elementNumber, 0) * f);
        }
    }
    return 2 * derivative;
}

double PadeJastrow::computeFirstEnergyDerivative() {
    m_particles          = m_system->getParticles();
    m_parameters         = m_system->getWeights();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += m_particles(i,j) * m_particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }
    return 0.5 * double(r.transpose()*r);
}

double PadeJastrow::computeSecondEnergyDerivative() {
    return 0;
}
