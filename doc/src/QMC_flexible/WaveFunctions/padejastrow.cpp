#include "padejastrow.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"

PadeJastrow::PadeJastrow(System* system, double alpha, double beta, Eigen::MatrixXd Gamma) :
        WaveFunction(m_system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 2;
    m_parameters.reserve(2);
    m_parameters.push_back(alpha);
    m_parameters.push_back(beta);
    m_Gamma = Gamma;
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
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);

    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }
    double PadeJastrowFactor = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        for(int j=0; j<i; j++) {
            R(i,j) = fabs(r(i) - r(j));
            PadeJastrowFactor += m_Gamma(i,j) * R(i,j)/(1 + m_parameters.at(1) * R(i,j));
        }
    }

    return exp(-0.5 * m_parameters.at(0) * (r.cwiseAbs2()).sum()) * PadeJastrowFactor;
}

double PadeJastrow::computeDerivative(Eigen::MatrixXd particles) {
    long m_numberOfParticles = particles.rows();
    long m_numberOfDimensions = particles.cols();
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m_numberOfParticles);
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(int i=0; i<m_numberOfParticles; i++) {
        double sqrtElementWise = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            sqrtElementWise += particles(i,j) * particles(i,j);
        }
        r(i) = sqrt(sqrtElementWise);
    }
    double derivative = -0.5 * m_parameters.at(0) * m_parameters.at(0) * (r.cwiseAbs2()).sum();
    for(int i=0; i<m_numberOfParticles; i++) {
        double derTerm = 0;
        for(int j=0; j<i; j++) {
            R(i,j) = fabs(r(i) - r(j));
            double f = 1/(1 + m_parameters.at(1) * R(i,j));
            derTerm += m_Gamma(i,j) * f * f;
            derivative += m_Gamma(i,j) * f * f * (m_parameters.at(1) * f + m_parameters.at(0) * r(i));
        }
        derivative -= 0.5 * derTerm * derTerm;
    }


    return derivative + 0.5 * m_parameters.at(0) * m_numberOfParticles * m_numberOfDimensions;
}

double PadeJastrow::computeEnergyDerivative(Eigen::MatrixXd particles) {
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
