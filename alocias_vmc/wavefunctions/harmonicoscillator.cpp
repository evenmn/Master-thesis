#ifdef HARMONICOSCILLATOR

#include "harmonicoscillator.h"
#include "../methods.h"
#include "../hermite/hermite.h"
#include "../slater.h"

#include <iostream>

Quantumdot::Quantumdot(Slater* sIn) : QuantumdotBasis::QuantumdotBasis() {
    /* set Slater(parent) object and switch interaction on by default and push
     * functions for derivatives with respect to variational parameters to
     * vector */
    slater = sIn;
    
    m_interaction = true;
    isFull = false;

    // add derivatives with respect to variational parameters to list in given
    // order
    variationalDerivativeFunctionList .
        push_back(&Quantumdot::variationalDerivativeExpression);
} // end constructor

Quantumdot::~Quantumdot() {
} // end deconstructor

void Quantumdot::initializeParameters(const double& w) {
    /* set number of particles and allocate space for matrices */
    omega = w;
    omegaSq = w * w;
    
    // initialize basis (wrapper)
    QuantumdotBasis::setup(slater->m_numParticles, slater->m_dim);

    // make sure Slater is full-shell (make sure n is a magic number)
    checkIfFullShell();
} // end function initializeParameters 

void Quantumdot::setParameters(const Eigen::VectorXd& newParameters) {
    /* update parameters */
    slater->m_parameters = newParameters;
    alpha = newParameters(0);
    aw = alpha*omega;
    sqaw = sqrt(aw);
} // end function setParameters

void Quantumdot::setInteraction(bool t) {
    /* toggle interaction */
    m_interaction = t;
} // end function setInteraction

void Quantumdot::checkIfFullShell() {
    /* make sure shell is filled */
    isFull = false;
    for (unsigned int i = 0; i < QuantumdotBasis::getMagic().size(); ++i) {
        if (slater->m_numParticles==QuantumdotBasis::getMagic(i)) {
            isFull = true;
            break;
        } // end if
    } // end fori
} // end function checkIfFullShell

double Quantumdot::variationalDerivativeExpression(const unsigned int& i, const
        unsigned int& j, const unsigned int& l) {
    /* calculate and return expression in derivative with respect to
     * variational parameter (derivative with respect to wavefunction in index
     * (i,j) in Slater matrix) */
    double sum = - 0.5 * omega * slater->getNewPosition(i).squaredNorm();
    for (unsigned int d = 0; d < slater->m_dim; ++d) {
        const int& n = QuantumdotBasis::getn(j,d);
        if (n==0) {
            continue;
        } // end if
        sum += n/alpha * m_SnewPositions(i,d) * m_hermite3DMatrix(i,d)(n-1) /
            m_hermite3DMatrix(i,d)(n);
    } // end ford
    return sum * slater->getWavefunction(i,j);
} // end function variationalDerivativeExpression

void Quantumdot::setHermite3DMatrix(const unsigned int& p) {
    for (unsigned int d = 0; d < slater->m_dim; ++d) {
        for (unsigned int j = 0; j < QuantumdotBasis::Cartesian::getn().size();
                ++j) {
            m_hermite3DMatrix(p,d)(j) = H(m_SnewPositions(p,d), j);
        } // end ford
    } // end forj
} // end function setHermite3DMatrix

void Quantumdot::set(const Eigen::MatrixXd& newPositions) {
    /* override for set function */
    m_SnewPositions = sqaw * newPositions;
    m_SoldPositions = m_SnewPositions;

    for (unsigned int i = 0; i < slater->getNumberOfParticles(); ++i) {
        setHermite3DMatrix(i);
    } // end fori

    m_oldHermite3DMatrix = m_hermite3DMatrix;
} // end function set;

std::string Quantumdot::setupDone() {
    /* return message of state of setup */
    std::string possibleN = " ";
    if (isFull) {
        /* return empty message if setup is successfull */
        return "";
    } else {
        /* return message if setup is successfull */
        for (unsigned int i = 0; i < QuantumdotBasis::getMagic().size(); ++i) {
            possibleN += std::to_string(QuantumdotBasis::getMagic(i)) + " ";
        } // end fori
        return "Slater not full, possible N:" + possibleN;
    } // end ifelse
} // end function setupDone

void Quantumdot::reSetAll() {
    /* function for reinitializing all matrices except alpha, beta, omega and
     * numParticles. Mainly used for testing */
    m_SnewPositions.setZero();
    m_SoldPositions.setZero();
} // end function reSetAll

void Quantumdot::initializeMatrices() {
    /* initialize matrices with default 0 */
    m_SnewPositions = Eigen::MatrixXd::Zero(slater->m_numParticles,
            slater->m_dim);
    m_SoldPositions = Eigen::MatrixXd::Zero(slater->m_numParticles,
            slater->m_dim);

    m_laplacianSumVec = Eigen::VectorXd::Zero(slater->m_numParticles/2);

    m_hermite3DMatrix = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
                      Eigen::Dynamic>::Constant(slater->m_numParticles,
                              slater->m_dim, Eigen::VectorXd::Zero(
                                  QuantumdotBasis::Cartesian::getn() .
                                  size()));
    m_oldHermite3DMatrix = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
                         Eigen::Dynamic>::Constant(slater->m_numParticles,
                                 slater->m_dim, Eigen::VectorXd::Zero(
                                     QuantumdotBasis::Cartesian::getn() .
                                     size()));
} // end function setMatricesToZero

void Quantumdot::update(const Eigen::VectorXd& newPosition, const unsigned int&
        p) {
    /* update position for particle p and calculate new wavefunction, keep old
     * values */
    m_SoldPositions.row(p) = m_SnewPositions.row(p);
    m_SnewPositions.row(p) = sqaw * slater->getNewPosition(p);

    m_oldHermite3DMatrix.row(p) = m_hermite3DMatrix.row(p);
    setHermite3DMatrix(p);
} // end function update

void Quantumdot::reset(const unsigned int& p) {
    /* set new position and wavefunction to old values for particle p */
    m_SnewPositions.row(p) = m_SoldPositions.row(p);

    m_hermite3DMatrix.row(p) = m_oldHermite3DMatrix.row(p);
} // end function reset

void Quantumdot::acceptState(const unsigned int& p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    m_SoldPositions.row(p) =  m_SnewPositions.row(p);

    m_oldHermite3DMatrix.row(p) = m_hermite3DMatrix.row(p);
} // end function acceptState

double Quantumdot::calculateWavefunction(const unsigned int& p, const unsigned
        int& j) {
    /* calculate and return new wavefunction for particle p in state j */
    double res = exp(-0.5*m_SnewPositions.row(p).squaredNorm());
    for (unsigned int d = 0; d < slater->m_dim; ++d) {
        res *= m_hermite3DMatrix(p,d)(QuantumdotBasis::getn(j,d));
    } // end ford
    return res;
} // end function calculateWavefunction

double Quantumdot::gradientExpression(const unsigned int& p, const int& j,
        const unsigned int& d) {
    /* calculate gradient expression */
    const int& n = QuantumdotBasis::getn(j,d);
    if (n==0) {
        return - sqaw * m_SnewPositions(p,d) * slater->getWavefunction(p,j);
    } else {
        return sqaw * (2*n * m_hermite3DMatrix(p,d)(n-1) /
                m_hermite3DMatrix(p,d)(n) - m_SnewPositions(p,d)) *
            slater->getWavefunction(p,j);
    } // end if
} // end function calculateGradient

const Eigen::VectorXd& Quantumdot::laplacianExpression(const unsigned int& i,
        const unsigned int& halfSize) {
    /* calculate and return expression involved in the laplacian */
    m_laplacianSumVec.fill(m_SnewPositions.row(i).squaredNorm() -
            slater->m_dim);
    for (unsigned int j = 0; j < halfSize; ++j) {
        for (unsigned int d = 0; d < slater->m_dim; ++d) {
            const int& n = QuantumdotBasis::getn(j,d);
            if (n==0) {
                continue;
            } else if (n==1) {
                m_laplacianSumVec(j) -= 4*n/m_hermite3DMatrix(i,d)(n) *
                    m_SnewPositions(i,d)*m_hermite3DMatrix(i,d)(n-1);
            } else {
                m_laplacianSumVec(j) += 4*n/m_hermite3DMatrix(i,d)(n) *
                    ((n-1)*m_hermite3DMatrix(i,d)(n-2) -
                     m_SnewPositions(i,d)*m_hermite3DMatrix(i,d)(n-1));
            } // end ifeifelse
        } // end ford
        m_laplacianSumVec(j) *= slater->getWavefunction(i,j);
    } // end forj
    m_laplacianSumVec *= aw;
    return m_laplacianSumVec;
} // end function laplacianExpression

double Quantumdot::potentialEnergy() {
    /* calculate and return potential energy */
    double P = 0.5 * omegaSq *
        slater->getNewPositions().rowwise().squaredNorm().sum();
    if (m_interaction) {
        /* run with interaction */
        if (slater->m_dim == 1) {
            /* remove singularity in case of 1d */
            for (unsigned int i = 0; i < slater->m_numParticles; ++i) {
                for (unsigned int j = i+1; j < slater->m_numParticles; ++j) {
                    P += 1. / sqrt((slater->getNewPosition(i) -
                                slater->getNewPosition(j)).cwiseAbs2().sum()
                            + 0.25*0.25);
                } // end forj
            } // end fori
        } else {
            /* 2d and 3d versions */
            for (unsigned int i = 0; i < slater->m_numParticles; ++i) {
                for (unsigned int j = i+1; j < slater->m_numParticles; ++j) {
                    P += 1. / slater->getNewDistance(i,j);
                } // end forj
            } // end fori
        } // end ifelse
    } // end if
    return P;
} // end function potentialEnergy

#endif
