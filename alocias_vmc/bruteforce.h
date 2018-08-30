#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H

#include "vmc.h"

template<class T>
class Bruteforce : public VMC<T> {
    private:
        void initializeSample() {
            /* calculate and set new positions, distances, wavefunction,
             * derivatives and parameters */
            VMC<T>::m_accumulativeValues = 0;

            Eigen::MatrixXd newPositions =
                Eigen::MatrixXd::Zero(VMC<T>::numParticles, VMC<T>::m_dim);
            for (unsigned int i = 0; i < VMC<T>::numParticles; ++i) {
                for (unsigned int j = 0; j < VMC<T>::m_dim; ++j) {
                    newPositions(i,j) = VMC<T>::stepMonteCarlo *
                        (VMC<T>::realDistribution(VMC<T>::genMT64) - 0.5);
                } // end forj
            } // end fori

            // set positions in wavefunction
            VMC<T>::wf->set(newPositions);
        } // end function initializeSample

        bool sample() {
            /* run Metropolis sampling, return true in case state is accepted
             * and false if rejected  */

            // move only 1 particle at a time
            unsigned int pIdx = VMC<T>::randomParticle();

            // find new positions
            Eigen::VectorXd newPosition(VMC<T>::m_dim);
            for (unsigned int d = 0; d < VMC<T>::m_dim; ++d) {
                newPosition(d) = VMC<T>::stepMonteCarlo *
                    (VMC<T>::realDistribution(VMC<T>::genMT64) - 0.5);
            } // end ford

            // set positions and wavefunction
            VMC<T>::wf->update(newPosition, pIdx);

            // calculate ratio (squared)
            double wavefunctionRatio = VMC<T>::MetropolisRatio(pIdx);
            wavefunctionRatio *= wavefunctionRatio;

            // update values according to Metropolis test
            if (VMC<T>::MetropolisTest(wavefunctionRatio)) {
                /* accept move, increment acceptance */
                VMC<T>::m_accumulativeValues.acceptance++;
                VMC<T>::wf->acceptState(pIdx);
                return true;
            } else {
                /* reject move, reset values */
                VMC<T>::wf->reset(pIdx);
                return false;
            } // end ifelse
        } // end function sample

    public:
        Bruteforce(T* wavefunction, const double step, const Eigen::VectorXd&
                initialParameters, const unsigned int maxIterations, int
                myRank, int numprocs) : VMC<T>(wavefunction, step,
                    initialParameters, maxIterations, myRank, numprocs) {
            /* initialize VMC with Slater object */
        } // end constructor

        virtual ~Bruteforce () {
        } // end deconstructor
};

#endif /* BRUTEFORCE_H */
