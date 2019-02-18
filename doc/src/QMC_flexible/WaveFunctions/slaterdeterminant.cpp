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
    double h = 1;
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

double SlaterDeterminant::evaluate(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    Eigen::MatrixXd D_up = updateMatrix(positions.head(m_numberOfFreeDimensions/2), H);
    Eigen::MatrixXd D_dn = updateMatrix(positions.tail(m_numberOfFreeDimensions/2), H);

    return D_up.determinant()*D_dn.determinant();
}

double SlaterDeterminant::evaluateSqrd(Eigen::VectorXd positions, Eigen::VectorXd radialVector, Eigen::MatrixXd distanceMatrix) {
    Eigen::MatrixXd D_up = updateMatrix(positions.head(m_numberOfFreeDimensions/2), H);
    Eigen::MatrixXd D_dn = updateMatrix(positions.tail(m_numberOfFreeDimensions/2), H);

    return D_up.determinant()*D_dn.determinant()*D_up.determinant()*D_dn.determinant();
}

double SlaterDeterminant::computeFirstDerivative(const Eigen::VectorXd positions, int k) {

    Eigen::MatrixXd dA_up = dA_matrix(positions.head(m_numberOfFreeDimensions/2));
    Eigen::MatrixXd dA_dn = dA_matrix(positions.tail(m_numberOfFreeDimensions/2));

    Eigen::MatrixXd D_up = updateMatrix(positions.head(m_numberOfFreeDimensions/2), H);
    Eigen::MatrixXd D_dn = updateMatrix(positions.tail(m_numberOfFreeDimensions/2), H);

    double Sum = 0;
    for(int j=0; j<m_numberOfParticlesHalf; j++) {
        if(k<m_numberOfFreeDimensions/2) {
            Sum += dA_up(k,j) * D_up.inverse()(j, int(k/2));
        }
        else {
            Sum += dA_up(k-m_numberOfFreeDimensions/2,j) * D_up.inverse()(j, int((k-m_numberOfFreeDimensions/2)/2));
        }
    }
    return Sum;
}

double SlaterDeterminant::computeSecondDerivative() {;
    return 0;
}

Eigen::VectorXd SlaterDeterminant::computeFirstEnergyDerivative(int k) {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}

Eigen::VectorXd SlaterDeterminant::computeSecondEnergyDerivative() {
    return Eigen::VectorXd::Zero(m_maxNumberOfParametersPerElement);
}
