#include "expnqs.h"
#include "../slaterjastrow.h"

ExpNQS::ExpNQS(SlaterJastrow* SJIn) : PadeJastrow(SJIn), RBMJastrow(SJIn) {
    SJ = SJIn;
} // end constructor

ExpNQS::~ExpNQS() {
} // end deconstructor

void ExpNQS::initializeMatrices(const unsigned int& displ) {
    /* initialize matrices in PJ and RBMJ (add 1 to displ sendt to RBM to
     * account for beta parameter in PJ */
    m_jastrowGradient3DMatrix = EigenMatVecXd::Constant(SJ->m_numParticles,
            SJ->m_numParticles, Eigen::VectorXd::Zero(SJ->m_dim));
    m_oldJastrowGradient3DMatrix = EigenMatVecXd::Constant(SJ->m_numParticles,
            SJ->m_numParticles, Eigen::VectorXd::Zero(SJ->m_dim));
    m_jastrowGradientVector = Eigen::VectorXd::Zero(SJ->m_dim);

    RBMJastrow::initializeMatrices(displ);

    m_jastrowGradientVector = Eigen::VectorXd::Zero(SJ->m_dim);

    for (auto i : RBMJastrow::jastrowVariationalDerivativeFunctionList) {
        jastrowVariationalDerivativeFunctionList.push_back(i);
    } // end fori
} // end function initializeMatrices

double ExpNQS::jastrowWavefunctionRatioExpression(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return log of wavefunction ratio between new and old for
     * index i,j */
    return PadeJastrow::spinFactor(i,j) * (SJ->getNewDistance(i,j) -
            SJ->getOldDistance(i,j));
} // end function jastrowWavefunctionRatioExpression

double ExpNQS::jastrowWavefunctionRatio() {
    /* return ratio */
    double res = 0;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            res += jastrowWavefunctionRatioExpression(i,j);
        } // end forj
    } // end fori
    return RBMJastrow::jastrowWavefunctionRatio() * exp(res);
} // end function jastrowWavefunctionRatio

double ExpNQS::jastrowWavefunctionRatio(const unsigned int& p) {
    /* return ratio for row p */
    double res = 0;
    for (unsigned int j = 0; j < p; ++j) {
        res += jastrowWavefunctionRatioExpression(p,j);
    } // end forj
    for (unsigned int j = p+1; j < SJ->m_numParticles; ++j) {
        res += jastrowWavefunctionRatioExpression(p,j);
    } // end forj

    return RBMJastrow::jastrowWavefunctionRatio(p) * exp(res);
} // end function jastrowWavefunctionRatio

const Eigen::VectorXd& ExpNQS::jastrowGradientExpression(const unsigned
        int& i, const unsigned int& j) {
    /* calculate and return log of ratio between gradient and wavefunction for
     * indices i,j */
    m_jastrowGradientVector = PadeJastrow::spinFactor(i,j) *
        SJ->getNewDistanceVector(i,j) / SJ->getNewDistance(i,j);
    return m_jastrowGradientVector;
} // end function jastrowGradientExpression

double ExpNQS::jastrowLaplacian() {
    /* calculate Laplacian */
    static double dimm1 = SJ->m_dim - 1;
    double res = RBMJastrow::jastrowLaplacian();
    double resExp = 0.0;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            resExp += PadeJastrow::spinFactor(i,j) * dimm1 /
                SJ->getNewDistance(i,j);
        } // end forj
    } // end fori
    res += 2*resExp;

    for (unsigned int k = 0; k < SJ->getNumberOfParticles(); ++k) {
        /* return gradient for row p */
        m_jastrowGradientVector = m_jastrowGradient3DMatrix.row(k).sum();

        res += m_jastrowGradientVector.squaredNorm();
        res += 2 * m_jastrowGradientVector.dot(RBMJastrow::gradient(k));
    } // end fork

    return res;
} // end function jastrowLaplacian

void ExpNQS::calculateGradient() {
    /* call updaters in PJ and RBMJ */
    m_oldJastrowGradient3DMatrix = m_jastrowGradient3DMatrix;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        for (unsigned int j = i+1; j < SJ->m_numParticles; ++j) {
            m_jastrowGradient3DMatrix(i,j) = jastrowGradientExpression(i,j);
            m_jastrowGradient3DMatrix(j,i) = -m_jastrowGradient3DMatrix(i,j); 
        } // end forj
        for (unsigned int j = 0; j < i; ++j) {
            m_jastrowGradient3DMatrix(i,j) = jastrowGradientExpression(i,j);
            m_jastrowGradient3DMatrix(j,i) = -m_jastrowGradient3DMatrix(i,j); 
        } // end forj
    } // end fori
    m_oldJastrowGradient3DMatrix = m_jastrowGradient3DMatrix;

    RBMJastrow::calculateGradient();
} // end function calculateGradient

void ExpNQS::calculateGradient(const unsigned int& p) {
    /* call updaters in PJ and RBMJ for row p */
    Methods::symSwapNoDiag(m_oldJastrowGradient3DMatrix,
            m_jastrowGradient3DMatrix, p, SJ->getSpan());
    
    for (unsigned int j = p+1; j < SJ->m_numParticles; ++j) {
        m_jastrowGradient3DMatrix(p,j) = jastrowGradientExpression(p,j);
        m_jastrowGradient3DMatrix(j,p) = -m_jastrowGradient3DMatrix(p,j); 
    } // end forj
    
    for (unsigned int j = 0; j < p; ++j) {
        m_jastrowGradient3DMatrix(p,j) = jastrowGradientExpression(p,j);
        m_jastrowGradient3DMatrix(j,p) = -m_jastrowGradient3DMatrix(p,j); 
    } // end forj

    RBMJastrow::calculateGradient(p);
} // end function calculateGradient

const Eigen::VectorXd& ExpNQS::gradient(const unsigned int& p) {
    /* calculate and return gradient for row p */
    m_jastrowGradientVector = RBMJastrow::gradient(p) +
        m_jastrowGradient3DMatrix.row(p).sum();

    return m_jastrowGradientVector;
} // end function gradient

void ExpNQS::resetGradient(const unsigned int& p) {
    /* revert to old gradient */
    Methods::symSwapNoDiag(m_oldJastrowGradient3DMatrix,
            m_jastrowGradient3DMatrix, p, SJ->getSpan());
} // end function resetGradient

void ExpNQS::acceptGradient(const unsigned int& p) {
    /* accept gradient */
    Methods::symSwapNoDiag(m_oldJastrowGradient3DMatrix,
            m_jastrowGradient3DMatrix, p, SJ->getSpan());
} // end function acceptGradient
