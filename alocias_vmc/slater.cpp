#include "slater.h"
#include "methods.h"

Slater::Slater(const unsigned int& dim, const unsigned int& n, const
        Eigen::VectorXd& initialParameters) : WF(this) {
    /* initialize dimensions used and number of particles n */
    m_dim = dim;
    m_numParticles = n;
    halfSize = m_numParticles / 2;

    // allocate space for initial parameters and derivatives (these are assumed
    // to include the Jastrow parts if given)
    WF::setParameters(initialParameters);
    m_firstDerivativesParameters =
        Eigen::VectorXd::Zero(initialParameters.size());
} // end function constructor

Slater::~Slater() {
} // end function deconstructor

const unsigned int& Slater::getSpan() const {
    /* return spanIdx */
    return spanIdx;
} // end function getSpan

const unsigned int& Slater::getDimension() const {
    /* return number of dimensions */
    return m_dim;
} // end function getDimension

const unsigned int& Slater::getNumberOfParticles() const {
    /* return number of rows in matrices */
    return m_numParticles;
} // end function getSize

const Eigen::VectorXd &Slater::getParameters() const {
    /* return vector containing derivatives with respect to parameters */
    return m_parameters;
} // end function getParameters

const Eigen::VectorXd &Slater::getVariationalDerivatives() const {
    return m_firstDerivativesParameters;
} // end functions set

const double &Slater::getNewPosition(const unsigned int& p, const unsigned
        int& d) const {
    /* return new position for particle p */
    return m_newPositions(p,d);
} // end function getNewPosition

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getNewPosition(const
        unsigned int& p) const {
    /* return new position for particle p */
    return m_newPositions.row(p);
} // end function getNewPosition

const Eigen::MatrixXd &Slater::getNewPositions() const {
    /* return matrix containing new positions */
    return m_newPositions;
} // end function getNewPosition

const double &Slater::getOldPosition(const unsigned int& p, const unsigned
        int& d) const {
    /* return new position for particle p */
    return m_oldPositions(p,d);
} // end function getOldPosition

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getOldPosition(const
        unsigned int& p) const {
    /* return old position for particle p */
    return m_oldPositions.row(p);
} // end function getOldPosition

const Eigen::MatrixXd &Slater::getOldPositions() const {
    /* return matrix containing old positions */
    return m_oldPositions;
} // end function getOldPosition

const double &Slater::getNewDistance(const unsigned int &p1, const unsigned
        int &p2) const {
    /* return new inter-electron distance for particle p1 to p2 */
    return m_newDistances(p1, p2);
} // end function getNewDistance

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getNewDistance(const
        unsigned int& p) const {
    /* return new inter-electron distance for particle p */
    return m_newDistances.row(p);
} // end function getNewDistance

const Eigen::MatrixXd &Slater::getNewDistanceMatrix() const {
    /* return matrix containing new inter-electron distances */
    return m_newDistances;
} // end function getNewDistance

const double &Slater::getOldDistance(const unsigned int &p1, const unsigned
        int &p2) const {
    /* return old inter-electron distance for particle p1 to p2 */
    return m_oldDistances(p1, p2);
} // end function getOldDistance

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getOldDistance(const
        unsigned int& p) const {
    /* return old inter-electron distance for particle p */
    return m_oldDistances.row(p);
} // end function getOldDistance

const Eigen::MatrixXd &Slater::getOldDistanceMatrix() const {
    /* return matrix containing old inter-electron distances */
    return m_oldDistances;
} // end function getOldDistanceMatrix

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getNewDistanceVector(const
        unsigned int& p1, const unsigned int& p2) {
    /* return element (p1,p2) in m_newDistanceMatrix */
    return m_newDistanceMatrix(p1,p2);
} // end function getNewDistanceVector

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getOldDistanceVector(const
        unsigned int& p1, const unsigned int& p2) {
    /* return element (p1,p2) in m_newDistanceMatrix */
    return m_oldDistanceMatrix(p1,p2);
} // end function getOldDistanceVector

const double &Slater::getGradient(const unsigned int& p, const unsigned
        int& d) const {
    /* return gradient for particle p in dimension d */
    return m_gradient(p,d);
} // end function getGradient

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getGradient(const
        unsigned int& p) const {
    /* return gradient for particle p */
    return m_gradient.row(p);
} // end function getGradient

const Eigen::MatrixXd &Slater::getGradientMatrix() const {
    /* return matrix containing gradient for all particles */
    return m_gradient;
} // end function getGradientMatrix

const Eigen::MatrixXd &Slater::getWavefunctionMatrix() const {
    /* return old Slater matrix */
    return m_newWavefunctionMatrix;
} // end function getOldWavefunctionMatrix

const Eigen::MatrixXd &Slater::getOldWavefunctionMatrix() const {
    /* return old Slater matrix */
    return m_oldWavefunctionMatrix;
} // end function getOldWavefunctionMatrix

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getWavefunction(const
        unsigned int& i) const {
    /* return Slater matrix element i,j */
    return m_newWavefunctionMatrix.row(i);
} // end function getWavefunction

const Eigen::Ref<const Eigen::RowVectorXd> Slater::getOldWavefunction(const
        unsigned int& i) const {
    /* return Slater matrix element i,j */
    return m_oldWavefunctionMatrix.row(i);
} // end function getWavefunction

const double &Slater::getWavefunction(const unsigned int& i, const unsigned
        int& j) const {
    /* return Slater matrix element i,j */
    return m_newWavefunctionMatrix(i,j);
} // end function getWavefunction

const double &Slater::getOldWavefunction(const unsigned int& i, const unsigned
        int& j) const {
    /* return old Slater matrix element i,j */
    return m_oldWavefunctionMatrix(i,j);
} // end function getOldWavefunction

const Eigen::MatrixXd &Slater::getInverseMatrix() const {
    /* return inverse Slater */
    return m_newWavefunctionMatrixInverse;
} // end function getNewInverseMatrix

const Eigen::MatrixXd &Slater::getOldInverseMatrix() const {
    /* return old inverse Slater */
    return m_oldWavefunctionMatrixInverse;
} // end function getOldInverseMatrix

const Eigen::Ref<const Eigen::MatrixXd> Slater::getInverse() const {
    /* return inverse Slater for current particle */
    return m_newWavefunctionMatrixInverse.block(waveIdx, 0, halfSize,
            halfSize);
} // end function getNewInverse

const Eigen::Ref<const Eigen::MatrixXd> Slater::getOldInverse() const {
    /* return old inverse Slater for current particle */
    return m_oldWavefunctionMatrixInverse.block(waveIdx, 0, halfSize,
            halfSize);
} // end function getOldInverse

double Slater::wavefunctionRatio() {
    /* return ratio between new- and old wavefunctions */
    return determinantRatio;
} // end function wavefunctionRatio

double Slater::wavefunctionRatio(const double& pIdx) {
    /* return ratio between new- and old wavefunctions (dummy for compliance
     * with general SlaterJastrow class */
    return determinantRatio;
} // end function wavefunctionRatio

void Slater::calculateWavefunction(const unsigned int& p) {
    /* calculate and set wavefunction for particle p */
    for (unsigned int j = 0; j < halfSize; ++j) {
        m_newWavefunctionMatrix(p,j) = WF::calculateWavefunction(p,j);
    } // end forj
} // end function calculateWavefunction

void Slater::setVariationalDerivatives() {
    /* set derivative with respect to variational parameter */
    if (variationalDerivativeFunctionList.size() != 0) {
        /* only calculate if function actually exists */
        for (unsigned int l = 0; l <
                variationalDerivativeFunctionList.size(); ++l) {
            m_firstDerivativesParameters(l) = 0;
            for (unsigned int i = 0; i < m_numParticles; ++i) {
                unsigned int jj, ii;
                if (i < halfSize) {
                    jj = 0;
                    ii = i;
                } else {
                    jj = halfSize;
                    ii = i - halfSize;
                } // end ifelse
                for (unsigned int j = 0; j < halfSize; ++j) {
                    m_firstDerivativesParameters(l) +=
                        (this->*(variationalDerivativeFunctionList .
                                 at(l)))(i,j,m_parameters(l))
                        * m_newWavefunctionMatrixInverse.block(jj, 0, halfSize,
                                halfSize)(j,ii);
                } // end forj
            } // end fori
        } // end forl
    } // end if
} // end function setVariationalDerivatives

void Slater::setWavefunction(const unsigned int& p) {
    /* set new wavefunction for particle p and keep old */
    m_oldWavefunctionMatrix.row(p) = m_newWavefunctionMatrix.row(p);
    m_oldWavefunctionMatrixInverse.block(waveIdx, 0, halfSize, halfSize) =
        m_newWavefunctionMatrixInverse.block(waveIdx, 0, halfSize, halfSize);

    // calculate and set new
    calculateWavefunction(p);

    // calculate ratio between new and old wavefunctions
    determinantRatio =
        Methods::determinantRatio(
                m_newWavefunctionMatrix.block(waveIdx, 0, halfSize, halfSize),
                m_oldWavefunctionMatrixInverse.block(waveIdx, 0, halfSize,
                    halfSize), invIdx);

    // calculate inverse(wavefunction)
    Methods::updateMatrixInverse<Eigen::MatrixXd>(
            m_oldWavefunctionMatrix.block(waveIdx, 0, halfSize, halfSize),
            m_newWavefunctionMatrix.block(waveIdx, 0, halfSize, halfSize),
            m_oldWavefunctionMatrixInverse.block(waveIdx, 0, halfSize,
                halfSize), m_newWavefunctionMatrixInverse.block(waveIdx, 0,
                    halfSize, halfSize), determinantRatio, invIdx);
} // end function setWavefunction

void Slater::setWavefunction() {
    /* set new wavefunction  and keep old */
    m_oldWavefunctionMatrix = m_newWavefunctionMatrix;
    m_oldWavefunctionMatrixInverse = m_newWavefunctionMatrixInverse;

    // calculate and set new
    for (unsigned int i = 0; i < m_numParticles; ++i) {
        for (unsigned int j = 0; j < halfSize; j++) {
            m_newWavefunctionMatrix(i,j) = WF::calculateWavefunction(i,j);
        } // end forj
    } // end fori

    // set ratio between new and old wavefunctions
    determinantRatio =
        (m_newWavefunctionMatrix.topRows(halfSize).determinant()
        * m_newWavefunctionMatrix.bottomRows(halfSize).determinant()) /
        (m_oldWavefunctionMatrix.topRows(halfSize).determinant()
         * m_oldWavefunctionMatrix.bottomRows(halfSize).determinant());

    // calculate inverse matrix
    m_newWavefunctionMatrixInverse.topRows(halfSize) =
        m_newWavefunctionMatrix.topRows(halfSize).inverse();
    m_newWavefunctionMatrixInverse.bottomRows(halfSize) =
        m_newWavefunctionMatrix.bottomRows(halfSize).inverse();
} // end function setWavefunction

void Slater::set(const Eigen::MatrixXd& newPositions) {
    /* set new positions and calculate new values for all matrices */
    reSetAll();

    m_newPositions = newPositions;
    m_oldPositions = m_newPositions;

    setDistances();
    m_oldDistances = m_newDistances;
    m_oldDistanceMatrix = m_newDistanceMatrix;

    wfset<WF, const Eigen::MatrixXd&>(static_cast<WF*>(this),newPositions);

    setWavefunction();
    m_oldWavefunctionMatrix = m_newWavefunctionMatrix;
    m_oldWavefunctionMatrixInverse = m_newWavefunctionMatrixInverse;

    determinantRatio = 1;
} // end function set

void Slater::setDistances() {
    /* set inter-electron distance matrix, keep old values */
    m_oldDistances = m_newDistances;
    m_oldDistanceMatrix = m_newDistanceMatrix;
    for (unsigned int i = 0; i < m_numParticles; ++i) {
        for (unsigned int j = i+1; j < m_numParticles; ++j) {
            m_newDistanceMatrix(i,j) = m_newPositions.row(i) -
                m_newPositions.row(j);
            m_newDistances(i,j) = m_newDistanceMatrix(i,j).norm();

            m_newDistanceMatrix(j,i) = - m_newDistanceMatrix(i,j);
            m_newDistances(j,i) = m_newDistances(i,j);
        } // end forj
    } // end fori
} // end function setDistances

void Slater::setDistances(const unsigned int &p) {
    /* set inter-electron distance matrix for particle p, keep old values */
    Methods::symSwapNoDiag(m_oldDistances, m_newDistances, p, spanIdx);
    Methods::symSwapNoDiag(m_oldDistanceMatrix, m_newDistanceMatrix, p,
            spanIdx);

    for (unsigned int j = p+1; j < m_numParticles; ++j) {
        m_newDistanceMatrix(p,j) = m_newPositions.row(p) -
            m_newPositions.row(j);
        m_newDistances(p,j) = m_newDistanceMatrix(p,j).norm();

        m_newDistanceMatrix(j,p) = - m_newDistanceMatrix(p,j);
        m_newDistances(j,p) = m_newDistances(p,j);
    } // end forj
    for (unsigned int j = 0; j < p; ++j) {
        m_newDistanceMatrix(p,j) = m_newPositions.row(p) -
            m_newPositions.row(j);
        m_newDistances(p,j) = m_newDistanceMatrix(p,j).norm();

        m_newDistanceMatrix(j,p) = - m_newDistanceMatrix(p,j);
        m_newDistances(j,p) = m_newDistances(p,j);
    } // end forj
} // end function setDistances

void Slater::setIndices(const unsigned int &p) {
    /* set indexes to sub-matrices in concatenated wavefunction- and inverse
     * matrices */
    if (p < halfSize) {
        invIdx = p;
        waveIdx = 0;
    } else {
        invIdx = p - halfSize;
        waveIdx = halfSize;
    } // end ifelse

    spanIdx = m_numParticles - (p+1);
} // end function setIndices

void Slater::setSize() {
    /* function for reinitializing all matrices except alpha, beta and omega.
     * Mainly used for testing */
    m_numParticles = getNumberOfParticles();
    halfSize = m_numParticles/2;
    initializeMatrices();
} // end function setSize 

void Slater::setSize(const unsigned int& n) {
    /* function for reinitializing all matrices except alpha, beta and omega.
     * Mainly used for testing */
    m_numParticles = n; 
    halfSize = m_numParticles/2;
    initializeMatrices();
} // end function setSize 

void Slater::reSetAll() {
    /* function for reinitializing all matrices except alpha, beta, omega and
     * numParticles. Mainly used for testing */
    m_newPositions.setZero();
    m_oldPositions.setZero();
    m_newDistances.setZero();
    m_oldDistances.setZero();
    Methods::matrix3DSetZero(m_oldDistanceMatrix);
    Methods::matrix3DSetZero(m_newDistanceMatrix);

    m_newWavefunctionMatrix.setZero();
    m_oldWavefunctionMatrix.setZero();
    m_newWavefunctionMatrixInverse.setZero();
    m_oldWavefunctionMatrixInverse.setZero();

    m_gradient.setZero();
    m_oldGradient.setZero();

    wfreSetAll<WF>(static_cast<WF*>(this));
} // end function reSetAll

void Slater::initializeMatrices() {
    /* initialize matrices with default 0 */
    m_newPositions = Eigen::MatrixXd::Zero(m_numParticles, m_dim);
    m_oldPositions = Eigen::MatrixXd::Zero(m_numParticles, m_dim);
    m_newDistances = Eigen::MatrixXd::Zero(m_numParticles, m_numParticles);
    m_oldDistances = Eigen::MatrixXd::Zero(m_numParticles, m_numParticles);

    m_oldDistanceMatrix = EigenMatVecXd::Constant(m_numParticles,
            m_numParticles, Eigen::VectorXd::Zero(m_dim));
    m_newDistanceMatrix = EigenMatVecXd::Constant(m_numParticles,
            m_numParticles, Eigen::VectorXd::Zero(m_dim));

    m_newWavefunctionMatrix = Eigen::MatrixXd::Zero(m_numParticles, halfSize);
    m_oldWavefunctionMatrix = Eigen::MatrixXd::Zero(m_numParticles, halfSize);
    m_newWavefunctionMatrixInverse = Eigen::MatrixXd::Zero(m_numParticles,
            halfSize);
    m_oldWavefunctionMatrixInverse = Eigen::MatrixXd::Zero(m_numParticles,
            halfSize);

    m_gradient = Eigen::MatrixXd::Zero(m_numParticles, m_dim);
    m_oldGradient = Eigen::MatrixXd::Zero(m_numParticles, m_dim);

    wfinitializeMatrices<WF>(static_cast<WF*>(this));
} // end function setMatricesToZero

void Slater::update(const Eigen::VectorXd& newPosition, const unsigned int&
        p) {
    /* update position for particle p and calculate new wavefunction, keep old
     * values */
    setIndices(p);

    m_oldPositions.row(p) = m_newPositions.row(p);
    m_newPositions.row(p) += newPosition;
    setDistances(p);

    wfupdate<WF, const Eigen::VectorXd&, const unsigned
        int&>(static_cast<WF*>(this),newPosition,p);

    setWavefunction(p);
} // end function update

void Slater::reset(const unsigned int& p) {
    /* set new position and wavefunction to old values for particle p */
    m_newPositions.row(p) = m_oldPositions.row(p);
    resetDistance(p);
    
    wfreset<WF, const unsigned int&>(static_cast<WF*>(this),p);

    m_newWavefunctionMatrix.row(p) = m_oldWavefunctionMatrix.row(p);
    m_newWavefunctionMatrixInverse.block(waveIdx, 0, halfSize, halfSize) =
        m_oldWavefunctionMatrixInverse.block(waveIdx, 0, halfSize, halfSize);

    determinantRatio = 1.;
} // end function reset

void Slater::resetDistance(const unsigned int& p) {
    Methods::symSwapNoDiag(m_newDistances, m_oldDistances, p, spanIdx);
    Methods::symSwapNoDiag(m_newDistanceMatrix, m_oldDistanceMatrix, p,
            spanIdx);
} // end function resetDistance

void Slater::resetGradient(const unsigned int& p) {
    /* set new gradient to old values for particle p */
    wfresetGradient<WF, const unsigned int&>(static_cast<WF*>(this),p);
    m_gradient.row(p) = m_oldGradient.row(p);
} // end function resetGradient

void Slater::acceptState(const unsigned int&p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    m_oldPositions.row(p) = m_newPositions.row(p);

    Methods::symSwapNoDiag(m_oldDistances, m_newDistances, p, spanIdx);
    Methods::symSwapNoDiag(m_oldDistanceMatrix, m_newDistanceMatrix, p,
            spanIdx);
    
    wfacceptState<WF, const unsigned int&>(static_cast<WF*>(this),p);

    m_oldWavefunctionMatrix.row(p) = m_newWavefunctionMatrix.row(p);
    m_oldWavefunctionMatrixInverse.block(waveIdx, 0, halfSize, halfSize) =
        m_newWavefunctionMatrixInverse.block(waveIdx, 0, halfSize, halfSize);
} // end function acceptState

void Slater::acceptGradient(const unsigned int&p) {
    /* accept state and set old gradient accordingly for particle
     * p */
    wfacceptGradient<WF, const unsigned int&>(static_cast<WF*>(this),p);
    m_oldGradient.row(p) = m_gradient.row(p);
} // end function acceptGradient

void Slater::calculateGradient() {
    /* calculate ratio between gradient of new wavefunction and new
     * wavefunction for all particles */
    m_oldGradient = m_gradient;
    m_gradient.setZero();
    for (unsigned int i = 0; i < m_numParticles; ++i) {
        unsigned int ii, wi;
        if (i < halfSize) {
            ii = i;
            wi = 0;
        } else {
            ii = i - halfSize;
            wi = halfSize;
        } // end ifelse
        for (unsigned int j = 0; j < halfSize; ++j) {
            for (unsigned int d = 0; d < m_dim; ++d) {
                m_gradient(i,d) += gradientExpression(i,j,d) *
                    m_newWavefunctionMatrixInverse.block(wi, 0, halfSize,
                            halfSize)(j,ii);
                } // end ford
            } // end forj
    } // end fori
    m_oldGradient = m_gradient;
} // end function gradient

void Slater::calculateGradient(const unsigned int& p) {
    /* calculate ratio between gradient of new wavefunction and new
     * wavefunction for particle p */
    m_oldGradient.row(p) = m_gradient.row(p);
    m_gradient.row(p).setZero();
    for (unsigned int d = 0; d < m_dim; ++d) {
        setGradient(p,d);
    } // end ford
} // end function gradient

void Slater::setGradient(const unsigned int& p, const unsigned int &d) {
    /* set gradient for particle p along dimension d */
    for (unsigned int j = 0; j < halfSize; ++j) {
        m_gradient(p,d) += gradientExpression(p,j,d) *
            m_newWavefunctionMatrixInverse.block(waveIdx, 0, halfSize,
                    halfSize)(j,invIdx);
    } // end forj
} // end function setGradient

double Slater::laplacian() {
    /* calculate ratio between laplacian of new wave function and new
     * wavefunction */
    double sum = 0.0;
    for (unsigned int i = 0; i < m_numParticles; ++i) {
        unsigned int ii, wi;
        if (i < halfSize) {
            ii = i;
            wi = 0;
        } else {
            ii = i - halfSize;
            wi = halfSize;
        } // end ifelse
        const Eigen::VectorXd& lapExpr = laplacianExpression(i, halfSize);
        for (unsigned int j = 0; j < halfSize; ++j) {
            sum += lapExpr(j) * m_newWavefunctionMatrixInverse.block(wi, 0,
                    halfSize, halfSize)(j,ii);
        } // end forj
    }  // end fori

    return sum;
} // end function laplacian

double Slater::kineticEnergy() {
    /* calculate and return kinetic energy */
    return 0.5 * laplacian();
} // end function kineticEnergy
