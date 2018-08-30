#ifndef IMPORTANCESAMPLING_H
#define IMPORTANCESAMPLING_H

#include "vmc.h"

class Slater;

template<class T>
class ImportanceSampling : public VMC<T> {
    private:
        Eigen::MatrixXd m_newQuantumForce, m_oldQuantumForce;
        
        std::normal_distribution<double> normalDistribution;

        void setQuantumForce() {
            /* set quantum force for all particles, keep old one */
            m_oldQuantumForce = m_newQuantumForce;
            for (unsigned int i = 0; i < VMC<T>::wf->getNumberOfParticles();
                    ++i) {
                m_newQuantumForce.row(i) = 2 * VMC<T>::wf->getGradient(i);
            } // end fori
        } // end function setQuantumForce

        void setQuantumForce(const unsigned int& p) {
            /* set quantum force for particle p, keep old one */
            m_oldQuantumForce.row(p) = m_newQuantumForce.row(p);
            m_newQuantumForce.row(p) = 2 * VMC<T>::wf->getGradient(p);
        } // end function setQuantumForce

        void setup() {
            /* allocate space for matrices when needed */
            if (!VMC<T>::m_hasSampleSetup) {
                /* allocate space only once */
                m_newQuantumForce = Eigen::MatrixXd::Zero(VMC<T>::numParticles,
                        VMC<T>::m_dim);
                m_oldQuantumForce = Eigen::MatrixXd::Zero(VMC<T>::numParticles,
                        VMC<T>::m_dim);

                VMC<T>::m_hasSampleSetup = true;
            } // end if
        } // end function setup

        void initializeSample() {
            /* set quantum force for all particles in case importance sampling is
             * switched on */
            setup();
            VMC<T>::m_accumulativeValues = 0;

            Eigen::MatrixXd newPositions =
                Eigen::MatrixXd::Zero(VMC<T>::numParticles, VMC<T>::m_dim);
            for (unsigned int i = 0; i < VMC<T>::numParticles; ++i) {
                for (unsigned int j = 0; j < VMC<T>::m_dim; ++j) {
                    newPositions(i,j) = normalDistribution(VMC<T>::genMT64) *
                        VMC<T>::sqrtStep;
                } // end forj
            } // end fori

            // set positions in wavefunction
            VMC<T>::wf->set(newPositions);

            // set gradient still if only Slater is used
            if constexpr (std::is_same<T, Slater>()) {
                /* figure it out at compile time */
                VMC<T>::wf->calculateGradient();
            } // end if

            // set quantum force
            m_oldQuantumForce.setZero();
            m_newQuantumForce.setZero();
            setQuantumForce();
            m_oldQuantumForce = m_newQuantumForce;
        } // end function initializeSample

        bool sample() {
            /* run Metropolis sampling, return true in case state is accepted and false
             * if rejected */

            // move only 1 particle at a time
            unsigned int pIdx = VMC<T>::randomParticle();

            // set positions and wavefunction
            Eigen::VectorXd newPosition(VMC<T>::m_dim);
            for (unsigned int d = 0; d < VMC<T>::m_dim; ++d) {
                newPosition(d) =
                    0.5*m_oldQuantumForce(pIdx,d)*VMC<T>::stepMonteCarlo +
                    normalDistribution(VMC<T>::genMT64) * VMC<T>::sqrtStep;
            } // end ford

            // Note: T=SlaterJastrow sets gradient as well
            VMC<T>::wf->update(newPosition, pIdx);

            // set gradient still if only Slater is used
            if constexpr (std::is_same<T, Slater>()) {
                /* figure it out at compile time */
                VMC<T>::wf->calculateGradient(pIdx);
            } // end if

            // calculate derivatives and set quantum force
            setQuantumForce(pIdx);

            // calculate ratio (squared)
            double wavefunctionRatio = VMC<T>::MetropolisRatio(pIdx);
            wavefunctionRatio *= wavefunctionRatio * greensFunction(pIdx);

            // update values according to Metropolis test
            if (VMC<T>::MetropolisTest(wavefunctionRatio)) {
                /* accept move, increment acceptance */
                VMC<T>::m_accumulativeValues.acceptance++;
                VMC<T>::wf->acceptState(pIdx);
                VMC<T>::wf->acceptGradient(pIdx);
                m_oldQuantumForce.row(pIdx) = m_newQuantumForce.row(pIdx);
                return true;
            } else {
                /* reject move, reset values */
                VMC<T>::wf->reset(pIdx);
                VMC<T>::wf->resetGradient(pIdx);
                m_newQuantumForce.row(pIdx) = m_oldQuantumForce.row(pIdx);
                return false;
            } // end ifelse
        } // end function sample

        double greensFunction(const unsigned int& pIdx) {
            /* calculate and return expression for greens function ratio used
             * in Metropolis test */
            return exp(0.25 * (m_oldQuantumForce.row(pIdx) +
                        m_newQuantumForce.row(pIdx)).dot(2 *
                        (VMC<T>::wf->getOldPosition(pIdx) -
                         VMC<T>::wf->getNewPosition(pIdx)) +
                        0.5*VMC<T>::stepMonteCarlo
                        * (m_oldQuantumForce.row(pIdx) -
                            m_newQuantumForce.row(pIdx)))) +
                (VMC<T>::numParticles - 1);
        } // end end function greensFunction
    
    public:
        ImportanceSampling(T* wavefunction, const double step, const
                Eigen::VectorXd& initialParameters, const unsigned int
                maxIterations, int myRank, int numprocs) : VMC<T>(wavefunction,
                    step, initialParameters, maxIterations, myRank, numprocs) {
            /* initialize normal distribution */
            normalDistribution = std::normal_distribution<double>{0,1};
        } // end constructor

        virtual ~ImportanceSampling () {
        } // end deconstructor
};

#endif /* IMPORTANCESAMPLING_H */
