#include "padejastrow.h"
#include "../slaterjastrow.h"
#include "../slater.h"

PadeJastrow::PadeJastrow(SlaterJastrow* sjin) {
    /* initialize dimnsions and number of particles */
    SJ = sjin;

    dimMinusOne = SJ->m_dim-1;

    // set spin factors that ensures cusp-conditions
    if (SJ->m_dim == 2) {
        /* 2D case */
        parallelSpinFactor = 1/3.;
        antiParallelSpinFactor = 1.0;
    } else if (SJ->m_dim == 3) {
        /* 3D case */
        parallelSpinFactor = 1/4.;
        antiParallelSpinFactor = 1/2.;
    } else {
        // FIXME: give error in case dim=1 or dim>3
    } // end ifeifelse
} // end constructor

PadeJastrow::~PadeJastrow() {
} // end deconstructor

void PadeJastrow::initializeMatrices(const unsigned int& displ) {
     /* initialize jastrow gradient vector and set displacement (by number of
      * variational parameters present in one-body Slater */
    m_jastrowGradient3DMatrix = EigenMatVecXd::Constant(SJ->m_numParticles,
            SJ->m_numParticles, Eigen::VectorXd::Zero(SJ->m_dim));
    m_oldJastrowGradient3DMatrix = EigenMatVecXd::Constant(SJ->m_numParticles,
            SJ->m_numParticles, Eigen::VectorXd::Zero(SJ->m_dim));

    m_jastrowGradientVector = Eigen::VectorXd::Zero(SJ->m_dim);

    m_parametersDispl = displ;

    jastrowVariationalDerivativeFunctionList .
        push_back(&PadeJastrow::jastrowVariationalDerivativeExpression);
} // end function initializeMatrices

double PadeJastrow::spinFactor(const unsigned int &i, const unsigned int &j) {
    /* return based on rule given in constructor (dependance on number of
     * dimension) */
    if (((i<SJ->halfSize) && (j<SJ->halfSize)) || ((i>=SJ->halfSize) &&
                (j>=SJ->halfSize))) {
        return parallelSpinFactor;
    } else {
        return antiParallelSpinFactor;
    } // end ifelse
} // end function spinFactor

double PadeJastrow::jastrowWavefunctionRatio() {
    /* calculate and return ratio between new- and old Jastrow factors */
    double res = 0;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            res += jastrowWavefunctionRatioExpression(i,j);
        } // end forj
    } // end fori
    return exp(res);
} // end function jastrowWavefunctionRatio

double PadeJastrow::jastrowWavefunctionRatio(const unsigned int& p) {
    /* calculate and return ratio between new- and old Jastrow factors */
    double res = 0;
    for (unsigned int j = 0; j < p; ++j) {
        res += jastrowWavefunctionRatioExpression(p,j);
    } // end forj
    for (unsigned int j = p+1; j < SJ->m_numParticles; ++j) {
        res += jastrowWavefunctionRatioExpression(p,j);
    } // end forj
    return exp(res);
} // end function jastrowWavefunctionRatio

double PadeJastrow::jastrowVariationalDerivativeExpression(const unsigned int&
        l) {
    /* calculate derivative with respect to beta */
    double res = 0.0;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            double denom = 1. / (SJ->m_parameters(m_parametersDispl) +
                    1./SJ->getNewDistance(i,j));
            res -= spinFactor(i,j) * denom*denom;
        } // end forj
    } // end fori
    return res;
} // end function jastrowVariationalDerivativeExpression

double PadeJastrow::jastrowWavefunctionRatioExpression(const unsigned int& i,
        const unsigned int& j) {
    /* calculate and return log of wavefunction ratio between new and old for
     * index i,j */
    return spinFactor(i,j) * (1./(1./SJ->getNewDistance(i,j) +
                SJ->m_parameters(m_parametersDispl)) -
            1./(1./SJ->getOldDistance(i,j) +
                SJ->m_parameters(m_parametersDispl)));
} // end function jastrowWavefunctionRatioExpression

const Eigen::VectorXd& PadeJastrow::jastrowGradientExpression(const unsigned
        int& i, const unsigned int& j) {
    /* calculate and return the ratio between gradient and wavefunction for
     * indices i,j */
    double denom = 1 + SJ->m_parameters(m_parametersDispl) *
        SJ->getNewDistance(i,j);
    m_jastrowGradientVector = SJ->getNewDistanceVector(i,j) * spinFactor(i,j) /
        (SJ->getNewDistance(i,j) * denom*denom);
    return m_jastrowGradientVector;
} // end function gradientExpression

double PadeJastrow::jastrowLaplacian() {
    /* calculate and return Laplacian of Jastrow factor */
    double res = 0;
    const double& beta = SJ->m_parameters(m_parametersDispl);
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            double denom = 1. / (1 + beta * SJ->getNewDistance(i,j));
            res += spinFactor(i,j) * (denom*denom) *
                (dimMinusOne/SJ->getNewDistance(i,j) - 2*beta * denom) +
                m_jastrowGradient3DMatrix(i,j).squaredNorm();
        } // end forj
    } // end fori
    return 2*res;
} // end function jastrowLaplacian

void PadeJastrow::calculateGradient() {
    /* calculate gradient and set 3D-matrix */
    m_oldJastrowGradient3DMatrix = m_jastrowGradient3DMatrix;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            m_jastrowGradient3DMatrix(i,j) = jastrowGradientExpression(i,j);
            m_jastrowGradient3DMatrix(j,i) = - m_jastrowGradient3DMatrix(i,j);
        } // end forj
        for (unsigned int j = 0; j < i; ++j) {
            m_jastrowGradient3DMatrix(i,j) = jastrowGradientExpression(i,j);
            m_jastrowGradient3DMatrix(j,i) = - m_jastrowGradient3DMatrix(i,j);
        } // end forj
    } // end fori
    m_oldJastrowGradient3DMatrix = m_jastrowGradient3DMatrix;
} // end function calculateGradient

void PadeJastrow::calculateGradient(const unsigned int& p) {
    /* calculate gradient for row p */
    Methods::symSwap(m_oldJastrowGradient3DMatrix, m_jastrowGradient3DMatrix,
            p);
    
    for (unsigned int j = p+1; j < SJ->m_numParticles; ++j) {
        m_jastrowGradient3DMatrix(p,j) = jastrowGradientExpression(p,j);
        m_jastrowGradient3DMatrix(j,p) = - m_jastrowGradient3DMatrix(p,j);
    } // end forj
    
    for (unsigned int j = 0; j < p; ++j) {
        m_jastrowGradient3DMatrix(p,j) = jastrowGradientExpression(p,j);
        m_jastrowGradient3DMatrix(j,p) = - m_jastrowGradient3DMatrix(p,j);
    } // end forj
} // end function calculateGradient

const Eigen::VectorXd& PadeJastrow::gradient(const unsigned int& p) {
    /* return gradient for row p */
    m_jastrowGradientVector = m_jastrowGradient3DMatrix.row(p).sum();

    return m_jastrowGradientVector;
} // end function gradient

void PadeJastrow::resetGradient(const unsigned int& p) {
    /* revert to old gradient */
//     Methods::symSwapNoDiag(m_oldJastrowGradient3DMatrix,
//             m_jastrowGradient3DMatrix, p, SJ->getSpan());
    Methods::symSwap(m_oldJastrowGradient3DMatrix, m_jastrowGradient3DMatrix,
            p);
} // end function resetGradient

void PadeJastrow::acceptGradient(const unsigned int& p) {
    /* accept gradient */
//     Methods::symSwapNoDiag(m_oldJastrowGradient3DMatrix,
//             m_jastrowGradient3DMatrix, p, SJ->getSpan());
    Methods::symSwap(m_oldJastrowGradient3DMatrix, m_jastrowGradient3DMatrix,
            p);
} // end function acceptGradient
