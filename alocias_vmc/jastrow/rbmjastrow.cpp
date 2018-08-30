#include "rbmjastrow.h"
#include "../slaterjastrow.h"
#include "../slater.h"
#include "../methods.h"

RBMJastrow::RBMJastrow(SlaterJastrow* sIn) {
    /* set parent object */
    SJ = sIn;
} // end constructor

RBMJastrow::~RBMJastrow() {
} // end deconstructor

void RBMJastrow::initializeMatrices(const unsigned int& displ) {
    /* initialize jastrow gradient vector and set displacement (by number of
     * variational parameters present in one-body Slater */
    m_parametersDispl = displ;
    m_numVisibleBias = SJ->m_numParticles*SJ->m_dim;
    m_numHiddenBias = (SJ->m_parameters.size() - displ - m_numVisibleBias) /
        (1+m_numVisibleBias);
    m_numWeights = m_numVisibleBias * m_numHiddenBias;
    
    m_hiddenDispl = m_parametersDispl + m_numVisibleBias;
    m_weightsDispl = m_hiddenDispl + m_numHiddenBias;

    m_weightsSumMatrix = Eigen::MatrixXd::Zero(SJ->m_numParticles,
            m_numHiddenBias);
    m_oldWeightsSumMatrix = Eigen::MatrixXd::Zero(SJ->m_numParticles,
            m_numHiddenBias);
    m_expSumVector = Eigen::VectorXd::Zero(m_numHiddenBias);
    m_oldExpSumVector = Eigen::VectorXd::Zero(m_numHiddenBias);

    // yes, this is messy af...
    for (unsigned int i = m_parametersDispl; i < m_hiddenDispl; ++i) {
            jastrowVariationalDerivativeFunctionList .
                push_back(&RBMJastrow::derivativeVisibleBias);
    } // end fori
    for (unsigned int i = m_hiddenDispl; i < m_weightsDispl; ++i) {
            jastrowVariationalDerivativeFunctionList .
                push_back(&RBMJastrow::derivativeHiddenBias);
    } // end fori
    for (unsigned int i = m_weightsDispl; i < SJ->m_parameters.size(); ++i) {
            jastrowVariationalDerivativeFunctionList .
                push_back(&RBMJastrow::derivativeWeight);
    } // end fori

    m_jastrowGradientVector = Eigen::VectorXd::Zero(SJ->m_dim);
} // end function initializeMatrices

unsigned int RBMJastrow::wi(const unsigned int& p,const unsigned int& j, const
        unsigned int& d) {
    /* return index for weight */
    return d + SJ->m_dim * (p + SJ->m_numParticles*j);
} // end function wi
        
void RBMJastrow::acceptGradient(const unsigned int& p) {
    /* set old sums and exponential sums to new one */
    m_oldWeightsSumMatrix.row(p) = m_weightsSumMatrix.row(p);
    m_oldExpSumVector = m_expSumVector;
} // end function acceptGradient

void RBMJastrow::resetGradient(const unsigned int& p) {
    /* set new sums and exponential sums to old one */
    m_weightsSumMatrix.row(p) = m_oldWeightsSumMatrix.row(p);
    m_expSumVector = m_oldExpSumVector;
} // end functon acceptGradient

void RBMJastrow::calculateGradient() {
    /* calculate sums and exponentials used in gradient */
    m_weightsSumMatrix.setZero();
    for (unsigned int p = 0; p < SJ->m_numParticles; ++p) {
        for (unsigned int j = 0; j < m_numHiddenBias; ++j) {
            for (unsigned int d = 0; d < SJ->m_dim; ++d) {
                m_weightsSumMatrix(p,j) += SJ->getNewPosition(p,d) *
                    SJ->m_parameters(m_weightsDispl + wi(p,j,d));
            } // end ford
        } // end forj
    } // end forp

    setWeightsSum();
    
    m_oldWeightsSumMatrix = m_weightsSumMatrix;
    m_oldExpSumVector = m_expSumVector;
} // end function calculateGradient

void RBMJastrow::calculateGradient(const unsigned int& p) {
    /* calculate sums and exponetial in row p */
    m_oldWeightsSumMatrix.row(p) = m_weightsSumMatrix.row(p);
    m_oldExpSumVector = m_expSumVector;

    for (unsigned int j = 0; j < m_numHiddenBias; ++j) {
        m_weightsSumMatrix(p,j) = 0;
        for (unsigned int d = 0; d < SJ->m_dim; ++d) {
            m_weightsSumMatrix(p,j) += SJ->getNewPosition(p,d) *
                SJ->m_parameters(m_weightsDispl + wi(p,j,d));
        } // end ford
    } // end forj

    setWeightsSum();
} // end function calculateGradient

void RBMJastrow::setWeightsSum() {
    /* set exponential with hidden bias and weights sum */
    for (unsigned int j = 0; j < m_numHiddenBias; ++j) {
        m_expSumVector(j) = exp(SJ->m_parameters(m_hiddenDispl + j) +
                m_weightsSumMatrix.col(j).sum());
    } // end forj
} // end function setWeightsSum

double RBMJastrow::jastrowWavefunctionRatio() {
    /* calculate ratio between new and old wavefunction */
    double oldBiasSum = 0.0;
    double newBiasSum = 0.0;
    for (unsigned int i = 0; i < SJ->m_numParticles; ++i) {
        const Eigen::VectorXd& visibleVec =
            SJ->m_parameters.segment(m_parametersDispl+i*SJ->m_dim, SJ->m_dim);
        newBiasSum += (SJ->getNewPosition(i).transpose() -
                visibleVec).squaredNorm();
        oldBiasSum += (SJ->getOldPosition(i).transpose() -
                visibleVec).squaredNorm();
    } // end fori
  
    // add 1 to each element in expSumVector (both of them) then take product
    return exp((oldBiasSum - newBiasSum)/2.) *
        (m_expSumVector.array()+1).prod() /
        (m_oldExpSumVector.array()+1).prod();
} // end function jastrowWavefunctionRatio

double RBMJastrow::jastrowWavefunctionRatio(const unsigned int& p) {
    /* calculate ratio between new and old wavefunction */
    const Eigen::VectorXd& visibleVec =
        SJ->m_parameters.segment(m_parametersDispl+p*SJ->m_dim, SJ->m_dim);
    double newBiasSum = (SJ->getNewPosition(p).transpose() -
            visibleVec).squaredNorm();
    double oldBiasSum = (SJ->getOldPosition(p).transpose() -
            visibleVec).squaredNorm();
   
    // add 1 to each element in expSumVector (both of them) then take product
    return exp((oldBiasSum - newBiasSum)/2.) *
        (m_expSumVector.array()+1).prod() /
        (m_oldExpSumVector.array()+1).prod();
} // end function jastrowWavefunctionRatio

double RBMJastrow::derivativeVisibleBias(const unsigned int& l) {
    /* calculate derivative with respect to visible bias l */
    unsigned int visibleIdx = l-m_parametersDispl;
    return SJ->getNewPosition(visibleIdx/SJ->m_dim, visibleIdx%(SJ->m_dim)) -
        SJ->m_parameters(l);
} // end function derivativeVisibleBias

double RBMJastrow::derivativeHiddenBias(const unsigned int& l) {
    /* calculate derivative with respect to hidden bias l */
    return 1. / (1. + 1./m_expSumVector(l-m_hiddenDispl));
} // end function derivativeHiddenBias

double RBMJastrow::derivativeWeight(const unsigned int& l) {
    /* calculate derivative with respect to weight l */
    unsigned int weightsIdx = l - m_weightsDispl;
    return SJ->getNewPosition((weightsIdx%(m_numVisibleBias))/SJ->m_dim,
            weightsIdx%SJ->m_dim) / (1 +
                1./m_expSumVector(weightsIdx/m_numVisibleBias));
} // end function derivativeWeight

const Eigen::VectorXd& RBMJastrow::gradient(const unsigned int& p) {
    /* calculate and return gradient for particle p */
    m_jastrowGradientVector =
        SJ->m_parameters.segment(m_parametersDispl+p*SJ->m_dim, SJ->m_dim) -
        SJ->getNewPosition(p).transpose();
    for (unsigned int j = 0; j < m_numHiddenBias; ++j) {
        double jfactor = 1. / (1 + 1./m_expSumVector(j));
        for (unsigned int d = 0; d < SJ->m_dim; ++d) {
            m_jastrowGradientVector(d) += SJ->m_parameters(m_weightsDispl +
                    wi(p,j,d)) * jfactor; 
        } // end ford
    } // end forj
    return m_jastrowGradientVector;
} // end function gradient

double RBMJastrow::jastrowLaplacian() {
    /* calculate and return Laplacian of Jastrow factor */
    double res = 0.0;
    for (unsigned int k = 0; k < SJ->m_numParticles; ++k) {
        m_jastrowGradientVector =
            SJ->m_parameters.segment(m_parametersDispl+k*SJ->m_dim, SJ->m_dim)
            - SJ->getNewPosition(k).transpose();
        for (unsigned int j = 0; j < m_numHiddenBias; ++j) {
            double expSum = m_expSumVector(j); 
            double expSum1 = 1. / (1 + expSum);
            double expSum1Neg = 1. / (1 + 1./expSum);
            for (unsigned int d = 0; d < SJ->m_dim; ++d) {
                const double& wkj = SJ->m_parameters(m_weightsDispl +
                        wi(k,j,d));
                res += wkj*wkj * expSum * expSum1*expSum1;
                m_jastrowGradientVector(d) += wkj * expSum1Neg;
            } // end ford
        } // end forj
        res += m_jastrowGradientVector.squaredNorm() - SJ->m_dim;
            
    } // end fork

    return res;
} // end function jastrowLaplacian
