#ifdef VARIATIONALHARTREEFOCKDOUBLEWELL

#include "variationalhartreefockdoublewell.h"
#include "../slater.h"
#include "../methods.h"
#include "../hermite/hermite.h"

#include <boost/math/special_functions/factorials.hpp>

VariationalHartreeFockDoubleWell::VariationalHartreeFockDoubleWell(Slater* sIn)
    :
        VariationalHartreeFockDoubleWellBasis::VariationalHartreeFockDoubleWellBasis()
{
    slater = sIn;
    m_interaction = true;
    
    // add derivatives with respect to variational parameters to list in given
    // order
    variationalDerivativeFunctionList .
        push_back(&VariationalHartreeFockDoubleWell::variationalDerivativeExpression);
} // end constructor

VariationalHartreeFockDoubleWell::~VariationalHartreeFockDoubleWell() {
} // end deconstructor

void VariationalHartreeFockDoubleWell::setHermite3DMatrix(const unsigned int&
        p) {
    for (unsigned int d = 0; d < slater->m_dim; ++d) {
        for (unsigned int j = 0; j <
                VariationalHartreeFockDoubleWellBasis::Cartesian::getn().
                size(); ++j) {
            m_hermite3DMatrix(p,d)(j) = H(m_SnewPositions(p,d), j);
        } // end ford
    } // end forj
} // end function setHermite3DMatrix

void VariationalHartreeFockDoubleWell::setParameters(const Eigen::VectorXd& newParameters) {
    /* update parameters */
    /* Set any extra parameters (i.e alpha) here if needed */
    slater->m_parameters = newParameters;

    aw = newParameters(0)*omega;
    sqrtaw = sqrt(aw);

    setHermiteNormalizations();
} // end function setParameters

void VariationalHartreeFockDoubleWell::setInteraction(bool t) {
    /* switch Coulomb interaction on/off */
    m_interaction = t;
} // end function setInteraction

void VariationalHartreeFockDoubleWell::initializeParameters(const double& w,
        Eigen::MatrixXd coefficientMatrix) {
    /* set number of particles and allocate space for matrices */
    omega = w;
    omegaSq = w * w;
    sqrtOmega = sqrt(w);

    // initialize basis (wrapper)
    VariationalHartreeFockDoubleWellBasis::setup(slater->m_dim);
    m_numBasis = VariationalHartreeFockDoubleWellBasis::DWC::C.rows();
    
    // set normalizations
    m_hermiteNormalizations =
        Eigen::ArrayXd::Zero(VariationalHartreeFockDoubleWellBasis::getn().
                size());

    m_C = coefficientMatrix.sparseView(1, 1e-6);
    m_C.makeCompressed();
    
    /* fill in any matrices dependant on sizes from basis here */

} // end function initializeParameters

double VariationalHartreeFockDoubleWell::variationalDerivativeExpression(const
        unsigned int& i, const unsigned int& j, const double& alphal) {
    /* calculate and return expression in derivative with respect to
     * variational parameter */
    /* fill in */
    double res = 0;
    double lenr = - 0.5 * omega * slater->getNewPosition(i).squaredNorm();
    for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt; ++lIt) {
        const double& lval = lIt.value();
        for (Eigen::SparseMatrix<double>::InnerIterator kIt(DWC::C, lIt.row());
                kIt; ++kIt) {
            double sum = lenr;
            for (unsigned int d = 0; d < slater->m_dim; ++d) {
                const int& n =
                    VariationalHartreeFockDoubleWellBasis::getn(kIt.row(),d);
                if (n==0) {
                    continue;
                } // end if
                sum += n/alphal * m_SnewPositions(i,d) *
                    m_hermite3DMatrix(i,d)(n-1) / m_hermite3DMatrix(i,d)(n);
            } // end ford
            res += lval * kIt.value() * sum *
                m_SWavefunctionMatrix(i,kIt.row());
        } // end forkIt
    } // end forlIt
    return res;
} // end function variationalDerivativeExpression

void VariationalHartreeFockDoubleWell::set(const Eigen::MatrixXd& newPositions) {
    /* set function */
    /* set any matrix/vector dependant on the position here. Is will be called
     * in Slater every time corresponding function there is called */
    m_SnewPositions = sqrtaw * newPositions;
    m_SoldPositions = m_SnewPositions;

    for (unsigned int i = 0; i < slater->getNumberOfParticles(); ++i) {
        setHermite3DMatrix(i);
    } // end fori
    m_oldHermite3DMatrix = m_hermite3DMatrix;

    setBasisWavefunction();
    m_SoldWavefunctionMatrix = m_SWavefunctionMatrix;
} // end function set

void VariationalHartreeFockDoubleWell::reSetAll() {
    /* function for reinitializing all matrices except alpha, m_parameters(1),
     * omega and numParticles. */
    /* reinitialize any matrix/vector dependant on the position here. Is will
     * be called in Slater every time corresponding function there is called */
    m_SnewPositions.setZero();
    m_SoldPositions.setZero();
    m_SWavefunctionMatrix.setZero();
    m_SoldWavefunctionMatrix.setZero();
} // end function reSetAll

void VariationalHartreeFockDoubleWell::initializeMatrices() {
    /* initialize matrices with default 0 */
    /* initialize any matrix/vector dependant on the position here. Is will be
     * called in Slater every time corresponding function there is called */
    m_laplacianSumVec = Eigen::VectorXd::Zero(slater->m_numParticles/2);
    
    m_SnewPositions = Eigen::MatrixXd::Zero(slater->m_numParticles,
            slater->m_dim);
    m_SoldPositions = Eigen::MatrixXd::Zero(slater->m_numParticles,
            slater->m_dim);
    m_SWavefunctionMatrix = Eigen::MatrixXd::Zero(slater->m_numParticles,
            m_numBasis);
    m_SoldWavefunctionMatrix = Eigen::MatrixXd::Zero(slater->m_numParticles,
            m_numBasis);

    m_hermite3DMatrix = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
                      Eigen::Dynamic>::Constant(slater->m_numParticles,
                              slater->m_dim, Eigen::VectorXd::Zero(
                                  VariationalHartreeFockDoubleWellBasis::
                                  Cartesian::getn().size()));
    m_oldHermite3DMatrix = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
                         Eigen::Dynamic>::Constant(slater->m_numParticles,
                                 slater->m_dim, Eigen::VectorXd::Zero(
                                     VariationalHartreeFockDoubleWellBasis::
                                     Cartesian::getn().size()));
} // end function setMatricesToZero

void VariationalHartreeFockDoubleWell::setHermiteNormalizations() {
    /* calcualate normalization factors for hermite functions */
    for (unsigned int i = 0; i < m_hermiteNormalizations.size(); ++i) {
        m_hermiteNormalizations(i) = sqrt(sqrtaw / (sqrt(M_PI) * pow(2,i)
                    * boost::math::factorial<double>(i)));
    } // end fori
} // end function setHermiteNormalizations

void VariationalHartreeFockDoubleWell::setBasisWavefunction() {
    /* set m_SWavefunctionMatrix */
    for (unsigned int i = 0; i < slater->m_numParticles; ++i) {
        setBasisWavefunction(i);
    } // end fori
} // end function setBasisWavefunction

void VariationalHartreeFockDoubleWell::setBasisWavefunction(const unsigned int&
        p) {
    /* set row p in m_SWavefunctionMatrix */
    double expFactor = exp(-0.5*m_SnewPositions.row(p).squaredNorm());
    for (unsigned int l = 0; l < m_numBasis; ++l) {
        m_SWavefunctionMatrix(p,l) = expFactor;
        for (unsigned int d = 0; d < slater->m_dim; ++d) {
            const int& n = VariationalHartreeFockDoubleWellBasis::getn(l,d);
            m_SWavefunctionMatrix(p,l) *= m_hermite3DMatrix(p,d)(n) *
                m_hermiteNormalizations(n);
        } // end ford
    } // end forl
} // end function setBasisWavefunction

void VariationalHartreeFockDoubleWell::update(const Eigen::VectorXd&
        newPosition, const unsigned int& p) {
    /* update position for particle p and calculate new wavefunction, keep old
     * values */
    /* update any matrix/vector dependant on the position or wavefunction here.
     * Is will be called in Slater every time corresponding function there is
     * called */
    m_SoldPositions.row(p) = m_SnewPositions.row(p);
    m_SnewPositions.row(p) = sqrtaw * slater->getNewPosition(p);
    
    m_oldHermite3DMatrix.row(p) = m_hermite3DMatrix.row(p);
    setHermite3DMatrix(p);

    m_SoldWavefunctionMatrix.row(p) = m_SWavefunctionMatrix.row(p);
    setBasisWavefunction(p);
} // end function update

void VariationalHartreeFockDoubleWell::reset(const unsigned int& p) {
    /* set new position and wavefunction to old values for particle p */
    /* reset(revert to old value) any matrix/vector dependant on the position
     * or wavefunction here.  Is will be called in Slater every time
     * corresponding function there is called */
    m_SnewPositions.row(p) = m_SoldPositions.row(p);
    m_SWavefunctionMatrix.row(p) = m_SoldWavefunctionMatrix.row(p);
    
    m_hermite3DMatrix.row(p) = m_oldHermite3DMatrix.row(p);

    m_SWavefunctionMatrix.row(p) = m_SoldWavefunctionMatrix.row(p);
} // end function reset

void VariationalHartreeFockDoubleWell::acceptState(const unsigned int&p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    /* accept(set old to current) any matrix/vector dependant on the the
     * wavefunction. Is will be called in Slater every time corresponding
     * function there is called */
} // end function acceptState

double VariationalHartreeFockDoubleWell::calculateWavefunction(const unsigned int& p, const unsigned
        int& j) {
    /* calculate and return new wavefunction for particle p in state j */
    /* fill in */
    double res = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt; ++lIt) {
        const double& lval = lIt.value();
        for (Eigen::SparseMatrix<double>::InnerIterator kIt(DWC::C, lIt.row());
                kIt; ++kIt) {
            res += lval * kIt.value() * m_SWavefunctionMatrix(p, kIt.row());
        } // end forkIt
    } // end forlIt
    return res;
} // end function calculateWavefunction

double VariationalHartreeFockDoubleWell::gradientExpression(const unsigned int& p, const int& j, const
        unsigned int& d) {
    /* calculate gradient expression */
    /* fill in */
    double res = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt; ++lIt) {
        const double& lval = lIt.value();
        for (Eigen::SparseMatrix<double>::InnerIterator kIt(DWC::C, lIt.row());
                kIt; ++kIt) {
            const int& nd =
                VariationalHartreeFockDoubleWellBasis::getn(kIt.row(), d);
            if (nd==0) {
                res -= lval * kIt.value() * m_SnewPositions(p,d) *
                    m_SWavefunctionMatrix(p, kIt.row());
            } else {
                res += lval * kIt.value() * (2*nd *
                        m_hermite3DMatrix(p,d)(nd-1) /
                        m_hermite3DMatrix(p,d)(nd) - m_SnewPositions(p,d)) *
                    m_SWavefunctionMatrix(p, kIt.row());
            } // end ifelse
        } // end forkIt
    } // end forlIt
    return res * sqrtaw;
} // end function calculateGradient

const Eigen::VectorXd&
VariationalHartreeFockDoubleWell::laplacianExpression(const unsigned int& i,
        const unsigned int& idx) {
    /* calculate and return expression involved in the laplacian */
    /* fill in */
    m_laplacianSumVec.setZero();
    double rlen = m_SnewPositions.row(i).squaredNorm() - slater->m_dim;
    for (unsigned int j = 0; j < idx; ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt;
                ++lIt) {
            const double& lval = lIt.value();
            for (Eigen::SparseMatrix<double>::InnerIterator kIt(DWC::C,
                        lIt.row()); kIt; ++kIt) {
                double lsum = rlen;
                for (unsigned int d = 0; d < slater->m_dim; ++d) {
                    const int& n =
                        VariationalHartreeFockDoubleWellBasis::getn(kIt.row(),d);
                    if (n==0) {
                        continue;
                    } else if (n==1) {
                        lsum -= 4*n/m_hermite3DMatrix(i,d)(n) *
                            m_SnewPositions(i,d)*m_hermite3DMatrix(i,d)(n-1);
                    } else {
                        lsum += 4*n/m_hermite3DMatrix(i,d)(n) *
                            ((n-1)*m_hermite3DMatrix(i,d)(n-2) -
                             m_SnewPositions(i,d)*m_hermite3DMatrix(i,d)(n-1));
                    } // end ifeifelse
                } // end ford
                m_laplacianSumVec(j) += lval * kIt.value() * lsum *
                    m_SWavefunctionMatrix(i, kIt.row());
            } // end forkIt
        } // end forlIt
    } // end forj
    m_laplacianSumVec *= aw;
    return m_laplacianSumVec;
} // end function laplacianExpression

double VariationalHartreeFockDoubleWell::potentialEnergy() {
    /* calculate and return potential energy */
    /* fill in */
    double P = 0.5 * omegaSq *
        slater->getNewPositions().rowwise().squaredNorm().sum() +
        1.0/8.0*omegaSq*4.0 * slater->m_numParticles;
    if (m_interaction) {
        /* run with interaction */
        for (unsigned int i = 0; i < slater->m_numParticles; ++i) {
            P -= 0.5*2.0*omegaSq * abs(slater->getNewPosition(i,0));
            for (unsigned int j = i+1; j < slater->m_numParticles; ++j) {
                P += 1. / slater->getNewDistance(i,j) ;
            } // end forj
        } // end fori
    } // end if
    return P;
} // end function potentialEnergy

#endif
