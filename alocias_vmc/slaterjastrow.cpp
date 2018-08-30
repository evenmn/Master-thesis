#include "slaterjastrow.h"
#include "methods.h"

SlaterJastrow::SlaterJastrow(const unsigned int& dim, const unsigned int& n,
        const Eigen::VectorXd& initialParameters) : Slater(dim, n,
            initialParameters), Jastrow(this) {
    /* initialize spin factors and matrices */

    m_numVariationalParameters = initialParameters.size();
} // end constructor

SlaterJastrow::~SlaterJastrow() {
} // end deconstructor
        
Eigen::RowVectorXd SlaterJastrow::getGradient(const unsigned int& p) {
    /* return the sum of gradient with respect to Jastrow factor and Slater
     * determinant for particle p */
    return Slater::m_gradient.row(p) + m_jastrowGradientMatrix.row(p);
} // end function getGradient

const double& SlaterJastrow::getJastrowGradient(const unsigned int& p, const
        unsigned int& d) const {
    /* return gradient for jastrow factor for particle p in dimension d */
    return m_jastrowGradientMatrix(p,d);
} // end function getJastrowGradient

const Eigen::Ref<const Eigen::RowVectorXd>
SlaterJastrow::getJastrowGradient(const unsigned int& p) const {
    /* return gradient for jastrow factor for particle p */
    return m_jastrowGradientMatrix.row(p);
} // end function getJastrowGradient

const Eigen::MatrixXd& SlaterJastrow::getJastrowGradientMatrix() const {
    /* return matrix containing gradient of jastrow factor for all particles */
    return m_jastrowGradientMatrix;
} // end function getJastrowGradientMatrix

const double &SlaterJastrow::getOldJastrowGradient(const unsigned int& p, const
        unsigned int& d) const {
    /* return gradient of jastrow factor for particle p in m_dimension d */
    return m_oldJastrowGradientMatrix(p,d);
} // end function getOldJastrowGradient

const Eigen::Ref<const Eigen::RowVectorXd>
SlaterJastrow::getOldJastrowGradient(const unsigned int& p) const {
    /* return old gradient of jastrow factor for particle p */
    return m_oldJastrowGradientMatrix.row(p);
} // end function getOldJastrowGradient

const Eigen::MatrixXd &SlaterJastrow::getOldJastrowGradientMatrix() const {
    /* return matrix containing old gradient of jastrow factor for all
     * particles */
    return m_oldJastrowGradientMatrix;
} // end function getOldJastrowGradientMatrix

void SlaterJastrow::setVariationalDerivatives() {
    /* calculate gradient with respect to variational parameters */
    for (unsigned int l = variationalDerivativeFunctionList.size(); l <
            m_numVariationalParameters; ++l) {
        m_firstDerivativesParameters(l) = (this ->*
                (jastrowVariationalDerivativeFunctionList .
                 at(l-variationalDerivativeFunctionList.size())))(l);
    } // end forl

    Slater::setVariationalDerivatives();
} // end function setVariationalDerivatives

double SlaterJastrow::wavefunctionRatio() {
    /* calculate and return ratio between new- and old Jastrow factors times
     * ratio of new and old Slater determinants */
    return jastrowWavefunctionRatio() * determinantRatio;
} // end function wavefunctionRatio

double SlaterJastrow::wavefunctionRatio(const unsigned int &p) {
    /* calculate and return ratio between new- and old Jastrow factors times
     * ratio of new and old Slater determinant */
    return jastrowWavefunctionRatio(p) * determinantRatio;
} // end function wavefunctionRatio

void SlaterJastrow::calculateGradient() {
    /* calculate gradient of Jastrow factor and Slater determinant */
    m_oldJastrowGradientMatrix = m_jastrowGradientMatrix;

    jastrowCalculateGradient<Jastrow>(static_cast<Jastrow*>(this));

    for (unsigned int i = 0; i < m_numParticles; ++i) {
        m_jastrowGradientMatrix.row(i) = Jastrow::gradient(i);
    } // end fori
    m_oldJastrowGradientMatrix = m_jastrowGradientMatrix;
    
    Slater::calculateGradient();
} // end function calculateGradient

void SlaterJastrow::calculateGradient(const unsigned int &p) {
    /* calculate gradient of Jastrow factor and Slater determinant for particle
     * p */
    m_oldJastrowGradientMatrix.row(p) = m_jastrowGradientMatrix.row(p);

    jastrowCalculateGradient<Jastrow,const unsigned
        int&>(static_cast<Jastrow*>(this), p);

    m_jastrowGradientMatrix.row(p) = Jastrow::gradient(p);

    Slater::calculateGradient(p);
} // end function calculateGradient

double SlaterJastrow::kineticEnergy() {
    /* return kinetic term (full expression) */
    return 0.5 * (Slater::laplacian() + jastrowLaplacian()) +
        Slater::m_gradient.cwiseProduct(m_jastrowGradientMatrix).sum();
} // end function kineticEnergy

void SlaterJastrow::initializeMatrices() {
    /* initialize matrices with default 0 */
    m_jastrowGradientMatrix = Eigen::MatrixXd::Zero(m_numParticles, m_dim);
    m_oldJastrowGradientMatrix = Eigen::MatrixXd::Zero(m_numParticles, m_dim);

    Jastrow::initializeMatrices(variationalDerivativeFunctionList.size());

    Slater::initializeMatrices();
} // end function initializeJastrowMatrices

void SlaterJastrow::set(const Eigen::MatrixXd newPositions) {
    /* set new positions and calculate new values for all matrices */
    Slater::set(newPositions);

    SlaterJastrow::calculateGradient();
} // end function set

void SlaterJastrow::acceptState(const unsigned int& p) {
    /* accept state and update wavefunction */
    SlaterJastrow::acceptGradient(p);

    Slater::acceptState(p);
} // end function acceptState

void SlaterJastrow::reset(const unsigned int& p) {
    /* reject state and revert to old wavefunction */
    SlaterJastrow::resetGradient(p);

    Slater::reset(p);
} // end function reset

void SlaterJastrow::update(const Eigen::VectorXd& newPosition, const unsigned
        int& p) {
    /* update state */
    Slater::update(newPosition, p);

    SlaterJastrow::calculateGradient(p);
} // end function update

void SlaterJastrow::resetGradient(const unsigned int& p) {
    /* set new jastrow gradient to old values for particle p */
    m_jastrowGradientMatrix.row(p) = m_oldJastrowGradientMatrix.row(p);

    jastrowResetGradient<Jastrow, const unsigned
        int&>(static_cast<Jastrow*>(this), p);

    Slater::resetGradient(p);
} // end function resetJastrowGradient

void SlaterJastrow::acceptGradient(const unsigned int& p) {
    /* accept state and set old gradient for jastrow and slater accordingly for
     * particle p */
    m_oldJastrowGradientMatrix.row(p) = m_jastrowGradientMatrix.row(p);

    jastrowAcceptGradient<Jastrow, const unsigned
        int&>(static_cast<Jastrow*>(this), p);

    Slater::acceptGradient(p);
} // end function acceptGradient
