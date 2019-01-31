#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

PadeJastrow::PadeJastrow(System* system, double beta, Eigen::MatrixXd Gamma) :
        WaveFunction(m_system) {
    assert(beta >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(beta);
    m_Gamma = Gamma;
    m_beta = beta;
}

double PadeJastrow::evaluate(Eigen::MatrixXd particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */

    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
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
            PadeJastrowFactor += m_Gamma(i,j) * R(i,j)/(1 + m_beta * R(i,j));
        }
    }

    return exp(PadeJastrowFactor);
}

double PadeJastrow::computeFirstDerivative(Eigen::MatrixXd particles, int k) {
    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
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

    double derivative = 0;
    for(int j=0; j<k; j++) {
        double f = 1/(1 + m_beta * R(k,j));
        derivative+= m_Gamma(k,j) * f * f;
    }
    return derivative;
}

double PadeJastrow::computeSecondDerivative(Eigen::MatrixXd particles) {
    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
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

    double derivative = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            double f = 1/(1 + m_beta * R(i,j));
            derivative += m_Gamma(i,j) * f * f * (1/r(i) - m_beta * f);
        }
    }
    return 2 * derivative;
}

double PadeJastrow::computeFirstEnergyDerivative(Eigen::MatrixXd particles) {
    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }
    return 0.5 * double(r.transpose()*r);
}

double PadeJastrow::computeSecondEnergyDerivative(Eigen::MatrixXd particles) {
    return 0;
}
