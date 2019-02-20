#include "slaterdeterminant.h"
#include <cassert>
#include <iostream>
#include "../system.h"

SlaterDeterminant::SlaterDeterminant(System* system,
                               int elementNumber) :
        WaveFunction(system) {
    m_elementNumber                     = elementNumber;
    m_numberOfFreeDimensions            = m_system->getNumberOfFreeDimensions();
    m_numberOfParticles                 = m_system->getNumberOfParticles();
    m_numberOfOrbitals                  = m_system->getNumberOfOrbitals();
    m_numberOfDimensions                = m_system->getNumberOfDimensions();
    m_numberOfParticlesHalf             = m_system->getNumberOfParticles()/2;
    m_maxNumberOfParametersPerElement   = m_system->getMaxNumberOfParametersPerElement();
    m_omega                             = m_system->getFrequency();
}

double H(double x, int n) {
    //Hermite polynomial of n'th degree

    if(n == 0) {
        return 1;
    }
    else if(n == 1) {
        return 2*x;
    }
    else {
        return 2*x*H(x,n-1)-2*(n-1)*H(x,n-2);
    }
}

double dH(double x, int n) {
    //Derivative of Hermite polynomial of n'th degree

    if(n == 0) {
        return 0;
    }
    else {
        return 2*n*H(x,n-1);
    }
}

void SlaterDeterminant::updateArrays(Eigen::VectorXd positions, int pRand) {
    m_oldPositions = m_positions;
    m_positions = positions;
    m_D_up = updateMatrix(m_positions.head(m_numberOfFreeDimensions/2), H);
    m_D_dn = updateMatrix(m_positions.tail(m_numberOfFreeDimensions/2), H);
    m_dD_up = dA_matrix(m_positions.head(m_numberOfFreeDimensions/2));
    m_dD_dn = dA_matrix(m_positions.tail(m_numberOfFreeDimensions/2));
}

void SlaterDeterminant::resetArrays() {
    m_positions = m_oldPositions;
}

void SlaterDeterminant::initializeArrays(Eigen::VectorXd positions) {
    m_positions = positions;
    m_D_up = updateMatrix(m_positions.head(m_numberOfFreeDimensions/2), H);
    m_D_dn = updateMatrix(m_positions.tail(m_numberOfFreeDimensions/2), H);
    m_dD_up = dA_matrix(m_positions.head(m_numberOfFreeDimensions/2));
    m_dD_dn = dA_matrix(m_positions.tail(m_numberOfFreeDimensions/2));
}

void SlaterDeterminant::updateParameters(Eigen::MatrixXd parameters) {
}

Eigen::MatrixXd SlaterDeterminant::list() {
    // Returns the index list used in Slater
    // For instance (0,0), (1,0), (0,1) for 6P in 2D
    //              (0,0,0), (1,0,0), (0,1,0), (0,0,1) for 6P in 3D etc..

    Eigen::MatrixXd order = Eigen::MatrixXd::Zero(m_numberOfParticlesHalf, m_numberOfDimensions);

    int counter = 0;
    // Two dimensions
    if (m_numberOfDimensions == 2) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            for(int s=i; s<m_numberOfOrbitals; s++) {
                int j = s - i;
                order(counter,1) = i;
                order(counter,0) = j;
                counter += 1;
            }
        }
    }

    // Three dimensions
    else if (m_numberOfDimensions == 3) {
        for(int i=0; i<m_numberOfOrbitals; i++) {
            for(int j=0; j<m_numberOfOrbitals; j++) {
                for(int s=i+j; s<m_numberOfOrbitals; s++) {
                    int k = s - i - j;
                    order(counter,0) = i;
                    order(counter,1) = j;
                    order(counter,2) = k;
                    counter += 1;
                }
            }
        }
    }
    return order;
}

void dA_element() {
    //Update element of dA
}


Eigen::VectorXd SlaterDeterminant::dA_row(const Eigen::VectorXd positions, int k) {
    //Update row of dA

    Eigen::VectorXd dA = Eigen::VectorXd::Zero(m_numberOfParticlesHalf);
    Eigen::MatrixXd order = list();

    // Find indices of relevant row
    Eigen::VectorXd a = Eigen::VectorXd::Zero(m_numberOfDimensions);
    int l = k%m_numberOfDimensions;
    for(int i=0; i<m_numberOfDimensions; i++) {
        a(i) = k-l+i;
    }

    // Find matrix
    for(int i=0; i<m_numberOfParticlesHalf; i++) {
        dA(i) = dH(positions(k), int(order(i, l)));
        for(int j=0; j<m_numberOfDimensions; j++) {
            if(int(a(j)) != k) {
                dA(i) *= H(positions(int(a(j))), int(order(i, j)));
            }
        }
    }
    return dA;
}


Eigen::MatrixXd SlaterDeterminant::dA_matrix(const Eigen::VectorXd positions) {
    //Initialize the entire dA matrix
    Eigen::MatrixXd dA = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticlesHalf);

    for(int k=0; k<m_numberOfParticlesHalf; k++) {
        dA.row(k) = dA_row(positions, k);
    }
    return dA;
}

double SlaterDeterminant::updateElement(Eigen::VectorXd positions, double basis(double, int), int i, int j) {
    // Updates an element in A-matrix

    Eigen::MatrixXd order = list();

    double element = 1;
    for(int k=0; k<m_numberOfDimensions; k++) {
        element *= basis(sqrt(m_omega) * positions(m_numberOfDimensions*i+k), int(order(j,k)));
    }

    return element;
}

Eigen::VectorXd SlaterDeterminant::updateRow(Eigen::VectorXd positions, double basis(double, int), int i) {
    // Updates a row in A-matrix

    Eigen::VectorXd A = Eigen::VectorXd::Ones(m_numberOfParticlesHalf);

    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        A(j) = updateElement(positions, basis, i, j);
    }
    return A;
}

Eigen::MatrixXd SlaterDeterminant::updateMatrix(Eigen::VectorXd positions, double basis(double, int)) {
    // Update the entire matrix

    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(m_numberOfParticlesHalf, m_numberOfParticlesHalf);

    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        A.row(j) = updateRow(positions, basis, j);
    }

    return A;
}

double SlaterDeterminant::evaluate() {
    return m_D_up.determinant() * m_D_dn.determinant();
}

double SlaterDeterminant::evaluateSqrd() {
    double WF = evaluate();
    return WF * WF;
}

double SlaterDeterminant::computeFirstDerivative(Eigen::VectorXd positions, int k) {
    Eigen::MatrixXd D_up = updateMatrix(positions.head(m_numberOfFreeDimensions/2), H);
    Eigen::MatrixXd D_dn = updateMatrix(positions.tail(m_numberOfFreeDimensions/2), H);
    Eigen::MatrixXd dD_up = dA_matrix(positions.head(m_numberOfFreeDimensions/2));
    Eigen::MatrixXd dD_dn = dA_matrix(positions.tail(m_numberOfFreeDimensions/2));
    double Sum = 0;
    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        if(k<m_numberOfFreeDimensions/2) {
            Sum += dD_up(k,j) * D_up.inverse()(j, int(k/2));
        }
        else {
            Sum += dD_dn(k-m_numberOfFreeDimensions/2,j) * D_dn.inverse()(j, int((k-m_numberOfFreeDimensions/2)/2));
        }
    }
    return Sum;
}

double SlaterDeterminant::computeSecondDerivative() {
    Eigen::VectorXd diff = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for(int k=0; k<m_numberOfFreeDimensions; k++) {
        for(int j=0; j<m_numberOfParticlesHalf; j++) {
            if(k<m_numberOfFreeDimensions/2) {
                diff(k) += m_dD_up(k,j) * m_D_up.inverse()(j, int(k/2));
            }
            else {
                diff(k) += m_dD_dn(k-m_numberOfFreeDimensions/2,j) * m_D_dn.inverse()(j, int((k-m_numberOfFreeDimensions/2)/2));
            }
        }
    }
    return -double(diff.transpose() * diff);
}

Eigen::VectorXd SlaterDeterminant::computeFirstEnergyDerivative(int k) {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}

Eigen::VectorXd SlaterDeterminant::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
