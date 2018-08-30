#include "padenqs.h"
#include "../slaterjastrow.h"

PadeNQS::PadeNQS(SlaterJastrow* SJIn) : PadeJastrow(SJIn), RBMJastrow(SJIn) {
    /* set parent object */
    SJ = SJIn;
} // end constructor

PadeNQS::~PadeNQS() {
} // end deconstructor

void PadeNQS::initializeMatrices(const unsigned int& displ) {
    /* initialize matrices in PJ and RBMJ (add 1 to displ sendt to RBM to
     * account for beta parameter in PJ */
    PadeJastrow::initializeMatrices(displ);
    RBMJastrow::initializeMatrices(displ + 1);

    m_jastrowGradientVector = Eigen::VectorXd::Zero(SJ->m_dim);

    for (auto i : PadeJastrow::jastrowVariationalDerivativeFunctionList) {
        jastrowVariationalDerivativeFunctionList.push_back(i);
    } // end fori
    for (auto i : RBMJastrow::jastrowVariationalDerivativeFunctionList) {
        jastrowVariationalDerivativeFunctionList.push_back(i);
    } // end fori
} // end function initializeMatrices

double PadeNQS::jastrowWavefunctionRatio() {
    /* return ratio */
    return PadeJastrow::jastrowWavefunctionRatio() *
        RBMJastrow::jastrowWavefunctionRatio();
} // end function jastrowWavefunctionRatio

double PadeNQS::jastrowWavefunctionRatio(const unsigned int& p) {
    /* return ratio for row  p */
    return PadeJastrow::jastrowWavefunctionRatio(p) *
        RBMJastrow::jastrowWavefunctionRatio(p);
} // end function jastrowWavefunctionRatio

double PadeNQS::jastrowLaplacian() {
    /* calculate Laplacian */
    double res = PadeJastrow::jastrowLaplacian() +
        RBMJastrow::jastrowLaplacian();
    for (unsigned int k = 0; k < SJ->getNumberOfParticles(); ++k) {
        res += 2 * PadeJastrow::gradient(k).dot(RBMJastrow::gradient(k));
    } // end fork
    return res;
} // end function jastrowLaplacian

void PadeNQS::calculateGradient() {
    /* call updaters in PJ and RBMJ */
    PadeJastrow::calculateGradient();
    RBMJastrow::calculateGradient();
} // end function calculateGradient

void PadeNQS::calculateGradient(const unsigned int& p) {
    /* call updaters in PJ and RBMJ for row p */
    PadeJastrow::calculateGradient(p);
    RBMJastrow::calculateGradient(p);
} // end function calculateGradient

const Eigen::VectorXd& PadeNQS::gradient(const unsigned int& p) {
    /* calculate and return gadient for row p */
    m_jastrowGradientVector = PadeJastrow::gradient(p) +
        RBMJastrow::gradient(p);
    return m_jastrowGradientVector;
} // end function gradient

void PadeNQS::resetGradient(const unsigned int& p) {
    /* revert to old gradient */
    PadeJastrow::resetGradient(p);
    RBMJastrow::resetGradient(p);
} // end function resetGradient

void PadeNQS::acceptGradient(const unsigned int& p) {
    /* accept gradient */
    PadeJastrow::acceptGradient(p);
    RBMJastrow::acceptGradient(p);
} // end function acceptGradient
