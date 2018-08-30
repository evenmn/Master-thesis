#include "hydrogenmolecule.h"
#include "../methods.h"
#include <iostream>

Hydrogenmolecule::Hydrogenmolecule(double protonDistance, const unsigned int
        dimension) {
    /* set distance between protons */
    dim = dimension;
    m_firstDerivativesParameters = Eigen::VectorXd(2);

    R = protonDistance;
    positionsProtons = Eigen::MatrixXd::Zero(2,dim);
    positionsProtons(0,0) = -protonDistance/2.;
    positionsProtons(1,0) = protonDistance/2.;
} // end constructor

Hydrogenmolecule::~Hydrogenmolecule() {
} // end deconstructor

const double &Hydrogenmolecule::getNewPosition(const unsigned int& p, const
        unsigned int& d) const {
    /* return new position for particle p */
    return m_newPositions(p,d);
} // end function getNewPosition

const Eigen::Ref<const Eigen::RowVectorXd>
Hydrogenmolecule::getNewPosition(const unsigned int& p) const {
    /* return new position for particle p */
    return m_newPositions.row(p);
} // end function getNewPosition

const Eigen::MatrixXd &Hydrogenmolecule::getNewPositions() const {
    /* return matrix containing new positions */
    return m_newPositions;
} // end function getNewPosition

const double &Hydrogenmolecule::getOldPosition(const unsigned int& p, const
        unsigned int& d) const {
    /* return old position for particle p */
    return m_oldPositions(p,d);
} // end function getOldPosition

const Eigen::Ref<const Eigen::RowVectorXd>
Hydrogenmolecule::getOldPosition(const unsigned int& p) const {
    /* return old position for particle p */
    return m_oldPositions.row(p);
} // end function getOldPosition

const Eigen::MatrixXd &Hydrogenmolecule::getOldPositions() const {
    /* return matrix containing old positions */
    return m_oldPositions;
} // end function getOldPosition

const unsigned int &Hydrogenmolecule::getNumberOfParticles() const {
    /* return number of particles (always 2 here) */
    return m_numParticles;
} // end function getNumberOfParticles

const double &Hydrogenmolecule::getGradient(const unsigned int& p, const
        unsigned int& d) const {
    /* return gradient for particle p in dimension d */
    return m_gradientMatrix(p,d);
} // end function getGradient

const Eigen::VectorXd Hydrogenmolecule::getGradient(const unsigned int& p)
    const {
    /* return gradient for particle p */
    return m_gradientMatrix.row(p);
} // end function getGradient

const Eigen::MatrixXd &Hydrogenmolecule::getGradientMatrix() const {
    /* return matrix containing gradient for both particles */
    return m_gradientMatrix;
} // end function getGradientMatrix

const double &Hydrogenmolecule::getJastrowGradient(const unsigned int& p, const
        unsigned int& d) const {
    /* return gradient with respect to Jastrow factor for particle p in
     * dimension d */
    return m_jastrowGradientMatrix(p,d);
} // end function getJastrowGradient

const Eigen::VectorXd Hydrogenmolecule::getJastrowGradient(const unsigned int&
        p) const {
    /* return gradient with respect to Jastrow factor for particle p */
    return m_jastrowGradientMatrix.row(p);
} // end function getJastrowGradient

const Eigen::MatrixXd &Hydrogenmolecule::getJastrowGradientMatrix() const {
    /* return matrix containing gradient with respect to Jastrow factor for
     * both particles */
    return m_jastrowGradientMatrix;
} // end function getJastrowGradientMatrix

const Eigen::VectorXd &Hydrogenmolecule::getVariationalDerivatives() const {
    /* return vector containing derivatives with respect to parameters */
    return m_firstDerivativesParameters;
} // end function getVariationalDerivatives

const unsigned int& Hydrogenmolecule::getDimension() const {
    return dim;
} // end function getDimension

void Hydrogenmolecule::setWavefunction(const unsigned int& p) {
    /* set new wavefunction for particle p and keep old */
    m_oldWavefunctionMatrix.row(p) = m_newWavefunctionMatrix.row(p);

    // calculate and set new
    for (unsigned int i = 0; i < 2; ++i) {
        m_newWavefunctionMatrix(p,i) = calculateWavefunction(p,i);
    } // end fori
} // end function setWavefunction

void Hydrogenmolecule::setWavefunction() {
    /* set new wavefunction and keep old */
    m_oldWavefunctionMatrix = m_newWavefunctionMatrix;

    // calculate and set new
    for (unsigned int p = 0; p < 2; ++p) {
        for (unsigned int i = 0; i < 2; ++i) {
            m_newWavefunctionMatrix(p,i) = calculateWavefunction(p,i);
        } // end fori
    } // end forp
} // end function setWavefunction

void Hydrogenmolecule::setParameters(const Eigen::VectorXd& newParameters) {
    /* update alpha and beta */
    alpha = newParameters(0);
    beta = newParameters(1);
} // end function setParameters

void Hydrogenmolecule::setNumParticles(const unsigned int n) {
    /* set number of particles (always 2 here) and allocate space for matrices
     * */
    m_numParticles = 2;

    m_newPositions = Eigen::MatrixXd::Zero(m_numParticles, dim);
    m_oldPositions = Eigen::MatrixXd::Zero(m_numParticles, dim);
    m_newDistances = Eigen::MatrixXd::Zero(m_numParticles, m_numParticles);
    m_oldDistances = Eigen::MatrixXd::Zero(m_numParticles, m_numParticles);
    m_newWavefunctionMatrix = Eigen::MatrixXd::Zero(m_numParticles,
            m_numParticles);
    m_oldWavefunctionMatrix = Eigen::MatrixXd::Zero(m_numParticles,
            m_numParticles);
    m_gradientMatrix = Eigen::MatrixXd::Zero(m_numParticles, dim);
    m_jastrowGradientMatrix = Eigen::MatrixXd::Zero(m_numParticles, dim);
} // end function setNumParticles

void Hydrogenmolecule::setVariationalDerivatives(bool jb) {
    /* calculate gradient with respect to parameters */

    // set derivative for alpha
    m_firstDerivativesParameters(0) = -
        ((m_newDistances(0,0)*m_newWavefunctionMatrix(0,0) -
          m_newDistances(0,1)*m_newWavefunctionMatrix(0,1)) /
         Methods::vecdiff(m_newWavefunctionMatrix.row(0)) -
         (m_newDistances(1,0)*m_newWavefunctionMatrix(1,0) -
          m_newDistances(1,1)*m_newWavefunctionMatrix(1,1)) /
         Methods::vecdiff(m_newWavefunctionMatrix.row(1)));

    // set derivative for beta
    if (jb) {
        m_firstDerivativesParameters(1) = -0.5 * pow(1/(beta +
                    1/(m_newPositions.row(0) - m_newPositions.row(1)).norm()),
                2);
    } // end if
} // end function setVariationalDerivatives

void Hydrogenmolecule::rset(const Eigen::MatrixXd& newPositions, unsigned int c)
{
    /* set new positions and calculate new distances and new wavefunction, call
     * again to set old matrices */
    m_oldPositions = m_newPositions;
    m_newPositions = newPositions;
    setDistances();
    setWavefunction();

    // break when new and old are set
    if (c < 1) {
        rset(newPositions, c+1);
    } else {
        return;
    } // end ifelse
} // end function set

void Hydrogenmolecule::set(const Eigen::MatrixXd& newPositions) {
    /* override for recursive set function */
    rset(newPositions, 0);
} // end function set

void Hydrogenmolecule::setDistances(const unsigned int& p) {
    /* set distance from particle p to protons (absolute norm) */
    for (unsigned int i = 0; i < 2; ++i) {
        m_newDistances(p,i) = (m_newPositions.row(p) -
                positionsProtons.row(i)).norm();
    } // end fori
} // end function setDistances

void Hydrogenmolecule::setDistances() {
    /* set distance from particles to protons */
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = 0; j < 2; ++j) {
            m_newDistances(i,j) = (m_newPositions.row(i) -
                    positionsProtons.row(j)).norm();
        } // end forj
    } // end fori
} // end function setDistances

void Hydrogenmolecule::update(const double& value, const unsigned int& p, const
        unsigned int& d) {
    /* update positions for particle p in dimension d and calculate new
     * distance and wavefunction, keep old values */
    m_oldPositions(p,d) = m_newPositions(p,d);
    m_newPositions(p,d) += value;
    setDistances(p);
    setWavefunction(p);
} // end function update

void Hydrogenmolecule::update(const Eigen::VectorXd& newPosition, const
        unsigned int& p) {
    /* update position for particle p and calculate new distance and
     * wavefunction, keep old values */
    m_oldPositions.row(p) = m_newPositions.row(p);
    m_newPositions.row(p) = newPosition;
    setDistances(p);
    setWavefunction(p);
} // end function update

void Hydrogenmolecule::reset(const unsigned int& p, const unsigned int& d) {
    /* set new positions, distances and wavefunction to old values for particle
     * p in dimension d*/
    m_newPositions(p,d) = m_oldPositions(p,d);
    setDistances(p);
    m_newWavefunctionMatrix.row(p) = m_oldWavefunctionMatrix.row(p);
} // end function reset

void Hydrogenmolecule::reset(const unsigned int& p) {
    /* set new positions, distances and wavefunction to old values for particle
     * p */
    m_newPositions.row(p) = m_oldPositions.row(p);
    setDistances(p);
    m_newWavefunctionMatrix.row(p) = m_oldWavefunctionMatrix.row(p);
} // end function reset

void Hydrogenmolecule::resetGradient(const unsigned int& p) {
    /* set new gradient to old values for particle p */
    gradient(p);
} // end function resetGradient

void Hydrogenmolecule::resetJastrowGradient(const unsigned int& p) {
    /* set new gradient to old values for particle p */
    calculateJastrowGradient(p);
} // end function resetGradient

void Hydrogenmolecule::acceptState(const unsigned int& p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    m_oldPositions.row(p) =  m_newPositions.row(p);
    m_oldWavefunctionMatrix.row(p) = m_newWavefunctionMatrix.row(p);
} // end function acceptState

void Hydrogenmolecule::acceptState(const unsigned int&p, const unsigned int& d)
{
    /* accept state and set old positions and matrices accordingly for particle
     * p along dimension d */
    m_oldPositions(p,d) =  m_newPositions(p,d);
    m_oldWavefunctionMatrix(p,d) = m_newWavefunctionMatrix(p,d);
} // end function acceptState

void Hydrogenmolecule::acceptGradient(const unsigned int& p) {
    /* accept state and set old gradient accordingly */
    gradient(p);
} // end function acceptGradient

void Hydrogenmolecule::acceptJastrowGradient(const unsigned int& p) {
    /* accept state and set old gradient accordingly */
    calculateJastrowGradient(p);
} // end function acceptGradient

double Hydrogenmolecule::calculateWavefunction(const unsigned int p) {
    /* calculate and return new wavefunction for particle p */
    return Methods::vecdiff((-alpha*m_newDistances.row(p)) .
            unaryExpr<double(*)(double)>(&exp));
} // end function calculateWavefunction

double Hydrogenmolecule::calculateWavefunction(const unsigned int p, const
        unsigned int pp) {
    /* calculate and return new wavefunction for particle p with proton pp */
    return exp(-alpha*m_newDistances(p,pp));
} // end function calculateWavefunction

double Hydrogenmolecule::wavefunctionRatio() {
    /* return ratio between new- and old wavefunctions */
    return Methods::rowwiseVecdiff(m_newWavefunctionMatrix).prod() /
        Methods::rowwiseVecdiff(m_oldWavefunctionMatrix).prod();
} // end function wavefunctionRatio

void Hydrogenmolecule::gradient() {
    /* calculate ratio between gradient of new wavefunction and new
     * wavefunction for particle p*/
    for (unsigned int i = 0; i < m_numParticles; ++i) {
        m_gradientMatrix.row(i) = alpha /
            Methods::vecdiff(m_newWavefunctionMatrix.row(i)) *
            ((m_newPositions.row(i) - positionsProtons.row(1)) *
             m_newWavefunctionMatrix(i,1) / m_newDistances(i,1) -
             (m_newPositions.row(i) - positionsProtons.row(0)) *
             m_newWavefunctionMatrix(i,0) / m_newDistances(i,0));
    }
} // end function gradient

void Hydrogenmolecule::gradient(const unsigned int& p) {
    /* calculate ratio between gradient of new wavefunction and new
     * wavefunction for particle p*/
    m_gradientMatrix.row(p) = alpha /
        Methods::vecdiff(m_newWavefunctionMatrix.row(p)) *
        ((m_newPositions.row(p) - positionsProtons.row(1)) *
         m_newWavefunctionMatrix(p,1) / m_newDistances(p,1) -
         (m_newPositions.row(p) - positionsProtons.row(0)) *
         m_newWavefunctionMatrix(p,0) / m_newDistances(p,0));
} // end function gradient

void Hydrogenmolecule::gradient(const unsigned int& p, const unsigned int& d) {
    /* calculate ratio between gradient of new wavefunction and new
     * wavefunction for particle p*/
    m_gradientMatrix.row(p) = alpha /
        Methods::vecdiff(m_newWavefunctionMatrix.row(p)) *
        ((m_newPositions.row(p) - positionsProtons.row(1)) *
         m_newWavefunctionMatrix(p,1) / m_newDistances(p,1) -
         (m_newPositions.row(p) - positionsProtons.row(0)) *
         m_newWavefunctionMatrix(p,0) / m_newDistances(p,0));
} // end function gradient

double Hydrogenmolecule::laplacian() {
    /* calculate ratio between laplacian of new wave function and new
     * wavefunction for particle p */
    double sum = 0;
    for (unsigned int i = 0; i < m_numParticles; ++i) {
        sum += alpha * ((dim + alpha) * (m_newWavefunctionMatrix(i,1) /
                    m_newDistances(i,1) - m_newWavefunctionMatrix(i,0) /
                    m_newDistances(i,0)) /
                Methods::vecdiff(m_newWavefunctionMatrix.row(i)) - alpha);
    } // end fori
    return sum;
} // end function laplacian

double Hydrogenmolecule::jastrow(const unsigned int &p) {
    /* calculate and return ratio between new- and old Jastrow factors */
    double r12New = (m_newPositions.row(0) - m_newPositions.row(1)).norm();
    double r12Old = (m_oldPositions.row(0) - m_oldPositions.row(1)).norm();
    return exp(0.5 * (r12New / (1 + beta*r12New) - r12Old / (1 +
                    beta*r12Old)));
} // end function jastrow

void Hydrogenmolecule::calculateJastrowGradient() {
    /* calculate gradient of Jastrow factor */
    Eigen::MatrixXd r12Vec = m_newPositions.row(0) - m_newPositions.row(1);
    double r12 = r12Vec.norm();
    double denom = 1 + beta*r12;
    m_jastrowGradientMatrix.row(0) = 0.5 * r12Vec / (r12 * denom * denom);
    m_jastrowGradientMatrix.row(1) = - m_jastrowGradientMatrix.row(0);
} // end function calculateJastrowGradient

void Hydrogenmolecule::calculateJastrowGradient(const unsigned int& p) {
    /* calculate gradient of Jastrow factor */
    calculateJastrowGradient();
} // end function calculateJastrowGradient

double Hydrogenmolecule::jastrowLaplacian() {
    /* calculate and return Laplacian of Jastrow factor */
    double r12 = (m_newPositions.row(0) - m_newPositions.row(1)).norm();
    double denom = 1 + beta*r12;
    return m_jastrowGradientMatrix.rowwise().squaredNorm().sum() + (dim - 1) /
        (r12 * denom * denom) - 2*beta / (denom * denom * denom);
} // end function jastrowLaplacian

double Hydrogenmolecule::potentialEnergy() {
    /* calculate and return potential energy */
    return -
        m_newDistances.unaryExpr(std::ptr_fun(Methods::inverseFraction<double,
                    double>)) . sum() + 1./(m_newPositions.row(0) -
                    m_newPositions.row(1)).norm() + 1./R;
//     return - 1./m_newDistances(0,0) - 1./m_newDistances(0,1) -
//         1./m_newDistances(1,0) - 1./m_newDistances(1,1) +
//         1./(m_newPositions.row(0) - m_newPositions.row(1)).norm() + 1./R;
} // end function potentialEnergy

double Hydrogenmolecule::kineticEnergy(bool j) {
    /* calculate and return kinetic energy */
    if (j) {
        return - (0.5 * (laplacian() + jastrowLaplacian()) +
                (m_gradientMatrix.cwiseProduct(m_jastrowGradientMatrix)) .
                sum());
    } else {
        return - 0.5 * laplacian();
    } // end ifelse
} // end function kineticEnergy
