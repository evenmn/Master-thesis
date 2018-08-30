#ifdef TESTWAVEFUNCTION

#include "testwavefunction.h"
#include "../slater.h"
#include "../methods.h"

TestWavefunction::TestWavefunction(Slater* sIn) : TestWavefunctionBasis::TestWavefunctionBasis() {
    slater = sIn;
    m_interaction = true;
} // end constructor

TestWavefunction::~TestWavefunction() {
} // end deconstructor

void TestWavefunction::setParameters(const Eigen::VectorXd& newParameters) {
    /* update parameters */
    /* Set any extra parameters (i.e alpha) here if needed */
    slater->m_parameters = newParameters;
} // end function setParameters

void TestWavefunction::setInteraction(bool t) {
    /* switch Coulomb interaction on/off */
    m_interaction = t;
} // end function setInteraction

void TestWavefunction::initializeParameters()
    /* set number of particles and allocate space for matrices */

    // initialize basis (wrapper)
    TestWavefunctionBasis::setup();
    
    /* fill in any matrices dependant on sizes from basis here */

} // end function initializeParameters

double TestWavefunction::calculateWavefunction(const unsigned int& p, const unsigned
        int& j) {
    /* calculate and return new wavefunction for particle p in state j */
    /* fill in */
    double res;
    return res;
} // end function calculateWavefunction

double TestWavefunction::gradientExpression(const unsigned int& p, const int& j, const
        unsigned int& d) {
    /* calculate gradient expression */
    /* fill in */
} // end function calculateGradient

const Eigen::VectorXd& TestWavefunction::laplacianExpression(const unsigned int& i, const
        unsigned int& idx) {
    /* calculate and return expression involved in the laplacian */
    /* fill in */
    m_laplacianSumVec.setZero();
    return m_laplacianSumVec;
} // end function laplacianExpression

double TestWavefunction::potentialEnergy() {
    /* calculate and return potential energy */
    /* fill in */
    double P = 0;
    return P;
} // end function potentialEnergy

double TestWavefunction::kineticEnergy() {
    /* calculate and return kinetic energy */
    return 0.5 * slater->laplacian();
} // end function kineticEnergy

#endif
