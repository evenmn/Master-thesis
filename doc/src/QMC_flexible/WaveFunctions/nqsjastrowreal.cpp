#include "nqsjastrowreal.h"
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include <iostream>

NQSJastrowReal::NQSJastrowReal(System* system, int elementNumber) :
        WaveFunction(system) {

    m_elementNumber                     = elementNumber;
    m_numberOfHiddenNodes               = m_system->getNumberOfHiddenNodes();
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    double sigma                        = m_system->getWidth();
    m_sigmaSqrd                         = sigma*sigma;
}

void NQSJastrowReal::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;
}

void NQSJastrowReal::resetArrays() {
    m_positions = m_oldPositions;
}

void NQSJastrowReal::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
}

void NQSJastrowReal::updateParameters(Eigen::MatrixXd parameters) {
    //m_a = (m_parameters.row(m_elementNumber)).head(m_numberOfFreeDimensions);
}

int fromWToParameterIndexR(int i, int j, int numberOfFreeDimensions) {
    return j*numberOfFreeDimensions + i;
}

//double fromParameterToWIndices(int i, int numberOfFreeDimensions) {
//    return i%numberPfFreeDimensions, int(i/numberOfFreeDimensions)
//}

Eigen::MatrixXd NQSJastrowReal::W() {
    m_parameters = m_system->getWeights();
    Eigen::VectorXd XXX = m_parameters.row(m_elementNumber).segment(m_numberOfHiddenNodes, m_numberOfFreeDimensions*m_numberOfHiddenNodes);
    Eigen::Map<Eigen::MatrixXd> W(XXX.data(), m_numberOfFreeDimensions, m_numberOfHiddenNodes);
    return W;
}

Eigen::VectorXd NQSJastrowReal::b() {
    m_parameters = m_system->getWeights();
    return (m_parameters.row(m_elementNumber)).head(m_numberOfHiddenNodes);
}

Eigen::VectorXd NQSJastrowReal::v(Eigen::VectorXd positions) {
    Eigen::VectorXd V = b() + W().transpose() * positions;
    return V;
}

Eigen::VectorXd NQSJastrowReal::f(Eigen::VectorXd positions) {
    Eigen::VectorXd V = v(positions);
    Eigen::VectorXd F = V.array().exp();

    //for(int i=0; i<m_numberOfHiddenNodes; i++) {
    //    F(i) = exp(V(i));
    //}

    return F;
}

Eigen::VectorXd NQSJastrowReal::g(Eigen::VectorXd positions) {
    Eigen::VectorXd F = f(positions);
    Eigen::VectorXd G = Eigen::VectorXd::Zero(m_numberOfHiddenNodes);

    G = (Eigen::VectorXd::Ones(m_numberOfHiddenNodes) + F).cwiseInverse();

    //for(int i=0; i<m_numberOfHiddenNodes; i++) {
    //    G(i) = 1/(1 + exp(V(i)));
    //}

    return G;
}

double NQSJastrowReal::evaluate() {
    //Eigen::VectorXd G = g(positions);
    //return (G.cwiseInverse()).prod();

    Eigen::VectorXd B = b();
    Eigen::MatrixXd w = W();
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        double Sum = 0;
        for(int i=0; i<m_numberOfFreeDimensions; i++) {
            Sum += w(i,j) * m_positions(i) / m_sigmaSqrd;
        }
        Prod *= 1 + exp(B(j) + Sum);
    }
    return Prod;
}

double NQSJastrowReal::evaluateSqrd() {
    //Eigen::VectorXd G = g(positions);
    //return ((G.cwiseInverse()).cwiseAbs2()).prod();

    Eigen::VectorXd B = b();
    Eigen::MatrixXd w = W();
    double Prod = 1;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        double Sum = 0;
        for(int i=0; i<m_numberOfFreeDimensions; i++) {
            Sum += w(i,j) * m_positions(i) / m_sigmaSqrd;
        }
        Prod *= 1 + exp(B(j) + Sum);
    }
    return Prod * Prod;
}

double NQSJastrowReal::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    //Eigen::VectorXd F = f(positions);
    //Eigen::VectorXd G = g(positions);
    //Eigen::MatrixXd w = W();
    //return double(w.row(k) * (F.cwiseProduct(G))) / m_sigmaSqrd;

    Eigen::VectorXd B = b();
    Eigen::MatrixXd w = W();
    double Sum1 = 0;
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        double Sum2 = 0;
        for(int i=0; i<m_numberOfFreeDimensions; i++) {
            Sum2 += w(i,j) * positions(i) / m_sigmaSqrd;
        }
        Sum1 += w(k,j) / (m_sigmaSqrd * (1 + exp(-B(j) - Sum2)));
    }
    return Sum1;
}

double NQSJastrowReal::computeSecondDerivative() {
    m_positions = m_system->getParticles();
    //Eigen::VectorXd F = f(m_positions);
    //Eigen::VectorXd G = g(m_positions);
    //Eigen::MatrixXd w = W();
    //return (w.cwiseAbs2() * F.cwiseProduct(G.cwiseAbs2())).sum();

    Eigen::VectorXd B = b();
    Eigen::MatrixXd w = W();
    double Sum1 = 0;
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            double Sum2 = 0;
            for(int i=0; i<m_numberOfFreeDimensions; i++) {
                Sum2 += w(i,j) * m_positions(i) / m_sigmaSqrd;
            }
            Sum1 += w(k,j) * w(k,j) / (m_sigmaSqrd * m_sigmaSqrd * (1 + exp(B(j) + Sum2)) * (1 + exp(-B(j) - Sum2)));
        }
    }
    return Sum1;
}

Eigen::VectorXd NQSJastrowReal::computeFirstEnergyDerivative(int k) {
    m_positions = m_system->getParticles();

    //Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);


    Eigen::VectorXd F = f(m_positions);
    Eigen::VectorXd G = g(m_positions);
    Eigen::MatrixXd w = W();
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    //gradients.head(m_numberOfHiddenNodes) = - 0.5 * (((w.row(k)).transpose()).cwiseProduct(G.cwiseProduct(F.cwiseAbs2()))) / m_sigmaSqrd;


    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int n = fromWToParameterIndexR(l, m, m_numberOfFreeDimensions);
            if(l == k) {
                gradients(n + m_numberOfHiddenNodes) = -0.5 * G(m) * F(m) * (1 + w(k,m) * G(m) * m_positions(k) / m_sigmaSqrd)/m_sigmaSqrd;
            }
            else {
                gradients(n + m_numberOfHiddenNodes) = -0.5 * w(k, m) * F(m) * G(m) * G(m) / (m_sigmaSqrd * m_sigmaSqrd);
            }
        }
    }


    Eigen::VectorXd B = b();
    //Eigen::MatrixXd w = W();

    // Update b
    for(int j=0; j<m_numberOfHiddenNodes; j++) {
        double Sum = 0;
        for(int i=0; i<m_numberOfFreeDimensions; i++) {
            Sum += w(i,j) * m_positions(i) / m_sigmaSqrd;
        }
        gradients(j) = - 0.5 * w(k,j) / (m_sigmaSqrd * (1 + exp(B(j) + Sum)) * (1 + exp(-B(j) - Sum)));
    }
    /*

    // Update W
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            double Sum = 0;
            for(int i=0; i<m_numberOfFreeDimensions; i++) {
                Sum += w(i,j) * m_positions(i) / m_sigmaSqrd;
            }
            int l = fromWToParameterIndex(k, j, m_numberOfFreeDimensions) + m_numberOfHiddenNodes;
            gradients(l) = -0.5 * w(k,j) * m_positions(k) * exp(B(j) + Sum) / (m_sigmaSqrd * m_sigmaSqrd * (1 + exp(B(j) + Sum) * (1 + exp(B(j) + Sum))));
        }
    }
    */
    return gradients;
}

Eigen::VectorXd NQSJastrowReal::computeSecondEnergyDerivative() {
    m_positions = m_system->getParticles();

    Eigen::VectorXd F = f(m_positions);
    Eigen::VectorXd G = g(m_positions);
    Eigen::MatrixXd w = W();
    Eigen::VectorXd B = b();
    Eigen::VectorXd gradients = Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);

    // Update b
    //gradients.head(m_numberOfHiddenNodes) = - 0.5 * ((w.cwiseAbs2().colwise().sum().transpose()).cwiseProduct(F.cwiseProduct((Eigen::VectorXd::Ones(m_numberOfHiddenNodes)-F).cwiseProduct(G.cwiseProduct(G.cwiseAbs2()))))) / (m_sigmaSqrd * m_sigmaSqrd);

    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfHiddenNodes; j++) {
            double Sum = 0;
            for(int i=0; i<m_numberOfFreeDimensions; i++) {
                Sum += w(i,j) * m_positions(i) / m_sigmaSqrd;
            }
            gradients(j) = - 0.5 * w(k,j) * w(k,j) * (1 - exp(B(j) - Sum)) / (m_sigmaSqrd * m_sigmaSqrd * (1 + exp(B(j) + Sum)) * (1 + exp(B(j) + Sum)) * (1 + exp(-B(j) - Sum)));
        }
    }

    // Update W
    for(int l=0; l<m_numberOfFreeDimensions; l++) {
        for(int m=0; m<m_numberOfHiddenNodes; m++) {
            int k = fromWToParameterIndexR(l, m, m_numberOfFreeDimensions);
            gradients(k + m_numberOfHiddenNodes) = F(m) * G(m) * G(m) * (2 * w(l,m) + w.cwiseAbs2().colwise().sum()(m) * m_positions(l) * (1 - F(m)) * G(m))/ (m_sigmaSqrd * m_sigmaSqrd);
        }
    }

    return gradients;
}
