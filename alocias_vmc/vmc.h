#ifndef VMC_H
#define VMC_H

#include <chrono>
#include <cmath>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "resampler.h"
#include "accumulate.h"

#include "methods.h"

template<class T> class Minimizer;

template<class T>
class VMC {
    friend class Minimizer<T>;
    private:
        unsigned int m_bins, thermalizationIdx;
        double m_rBinsMax, m_rStep;

        Eigen::ArrayXd m_histogram;
        Eigen::ArrayXXd m_radial2DHistogram;
        Eigen::Array<Eigen::ArrayXd, Eigen::Dynamic, Eigen::Dynamic>
            m_radial3DHistogram;

        Eigen::VectorXd m_oldParameters, m_newDerivativeParameters,
            m_oldDerivativeParameters;

        std::string m_oneBodyDensityFilename;
        std::string m_radialDensityFilename;

        std::uniform_int_distribution<unsigned int> discreteDistribution;

        double tmpPotentialEnergy, tmpKineticEnergy, tmpEnergy;

        bool m_findStep, m_makeOnebodyDensities, m_makeRadialDensities;

        void derivativeVariational() {
            /* calculate derivative of expectation value */
            m_newDerivativeParameters = 2. *
                (m_accumulativeValues.psiDerivativeTimesEnergy -
                 m_accumulativeValues.energy *
                 m_accumulativeValues.psiDerivative);
// 
//             m_newDerivativeParameters = 2./m_maxIterations *
//                 (m_accumulativeValues.psiDerivativeTimesEnergySquared +
//                  m_accumulativeValues.psiDerivativeTimesEnergy -
//                  (m_accumulativeValues.energySquared +
//                   m_accumulativeValues.energy)  *
//                  m_accumulativeValues.psiDerivative);
        } // end function derivativeEnergy

        void writeRadialHistogramToFile() {
            /* write histogram for radial densities to file */
            if (m_dim == 2) {
                std::ofstream outfile;
                outfile.open(m_radialDensityFilename);
                outfile << m_radial2DHistogram;
                outfile.close();
            } else if (m_dim == 3) {
                std::ofstream outfile;
                outfile.open(m_radialDensityFilename);
                outfile << m_radial3DHistogram;
                outfile.close();
            } // end ifeif
        } // end function write3DHistogramToFile

        void writeHistogramToFile() {
            /* write histogram for one-body densities to file */
            std::ofstream outfile;
            outfile.open(m_oneBodyDensityFilename);
            outfile << m_histogram;
            outfile.close();
        } // end function writeHistogramToFile

        void accumulateRadialDensities() {
            for (unsigned int p = 0; p < wf->getNumberOfParticles(); ++p) {
                for (unsigned int b = 0; b < m_bins; ++b) {
                    if ((fabs(wf->getNewPosition(p,0)) >= b*m_rStep) &&
                            (fabs(wf->getNewPosition(p,0)) < (b+1)*m_rStep)) {
                        for (unsigned int a = 0; a < m_bins; ++a) {
                            if ((fabs(wf->getNewPosition(p,1)) >= a*m_rStep) &&
                                    (fabs(wf->getNewPosition(p,1)) <
                                     (a+1)*m_rStep)) {
                                m_radial2DHistogram(a,b) +=1;
                            } // end if
                        } // end fora
                    } // end if
                } // end forb
            } // end forp
        } // and funtion acumulateRadialDensities

        void accumulateOnebodyDensities() {
            /* increment histogram values */
            for (unsigned int p = 0; p < wf->getNumberOfParticles(); ++p) {
                double rpNorm = wf->getNewPosition(p).norm();
                for (unsigned int b = 0; b < m_bins; ++b) {
                    if ((rpNorm >= b*m_rStep) && (rpNorm < (b+1)*m_rStep)) {
                        m_histogram(b) += 1;
                    } // end if
                } // end forb
            } // end forp
        } // end function accumulateOnebodyDensities

        void calculateLocalValues() {
            /* calculate local energy, kinetic and potential parts */
            tmpPotentialEnergy = wf->potentialEnergy();
            tmpKineticEnergy = -wf->kineticEnergy();
            tmpEnergy = tmpPotentialEnergy + tmpKineticEnergy;
        } // end function calculateLocalValues

        void accumulateLocalValues() {
            /* sum up values */
            m_accumulativeValues.energy += tmpEnergy;
            m_accumulativeValues.energySquared += tmpEnergy*tmpEnergy;
            m_accumulativeValues.potentialEnergy += tmpPotentialEnergy;
            m_accumulativeValues.kineticEnergy += tmpKineticEnergy;
        } // end function accumulateLocalValues

        void accumulateLocalDerivatives() {
            /* calculate and accumulate local values for derivatives */
            m_accumulativeValues.psiDerivative +=
                wf->getVariationalDerivatives();
            m_accumulativeValues.psiDerivativeTimesEnergy +=
                wf->getVariationalDerivatives() * tmpEnergy;
            m_accumulativeValues.psiDerivativeTimesEnergySquared +=
                wf->getVariationalDerivatives() * tmpEnergy * tmpEnergy;
        } // end function accumulateLocalDerivatives

        void accumulateSample(const unsigned int& i, const bool& accepted) {
            /* sum up values, only calculate new values if state is accepted
             * (that is if accepted is true) */
            
            if (i < thermalizationIdx) {
                return;
            } else if (i == thermalizationIdx) {
                m_accumulativeValues.acceptance = (accepted ? 1 :
                        0);
            } // end if

            if (accepted || (i==thermalizationIdx)) {
                /* calculate local energy and save kinetic and potential parts
                 * along with derivatives with respect to variational
                 * parameters */
                calculateLocalValues();
                wf->setVariationalDerivatives();

                accumulateLocalValues();
                accumulateLocalDerivatives();
            } else if (!accepted && (i>0)) {
                /* use old energies and derivatves in case of rejection */
                accumulateLocalValues();
                accumulateLocalDerivatives();
            } // end ifeifeifeif

            if (m_makeOnebodyDensities) {
                accumulateOnebodyDensities();
            } // end if

            if (m_makeRadialDensities) {
                accumulateRadialDensities();
            } // end if

            // gather values for resampling
            if (rs) {
                /* gather values for resampling if resampler object is given */
                rs->samples(i-thermalizationIdx) = tmpEnergy;
            } // end if
        } // end function accumulateSample

        int m_rank, m_numprocs;

        Resampler *rs;

    public:
        VMC(T* wavefunction, const double step, const Eigen::VectorXd&
                initialParameters, const unsigned int maxIterations, int
                myRank, int numprocs) {
            /* initialize Slater object */
            wf = wavefunction;

            thermalizationIdx = 50000;

            m_makeOnebodyDensities = false;
            m_makeRadialDensities = false;

            m_dim = wf->getDimension();
            m_findStep = false;
            m_hasSampleSetup = false;
            stepMonteCarlo = step;
            sqrtStep = sqrt(stepMonteCarlo);
            numParticles = wf->getNumberOfParticles();
           
            std::mt19937_64 genMT64(std::stol(std::to_string(
                            std::chrono::high_resolution_clock::now().
                            time_since_epoch().count()).substr(10)));
            realDistribution = std::uniform_real_distribution<double>{0, 1};
            discreteDistribution = std::uniform_int_distribution<unsigned int>{0,
                    numParticles-1};

            numParameters = initialParameters.size();
            wf->setParameters(initialParameters);
            m_accumulativeValues << numParameters;
            m_maxIterations = maxIterations;
            m_rank = myRank;
            m_numprocs = numprocs;

            rs = NULL;
        } // end constructor

        virtual ~VMC () {
        } // end deconstructor
        
        T* wf;

        double sampler(const Eigen::VectorXd& parameters, const unsigned int
                progressDivider=0) {
            /* overload which takes in and sets new parameters */
            wf->setParameters(parameters);
            return sampler(progressDivider);
        } // end function sampler

        double sampler(const unsigned int progressDivider=0) {
            /* run Monte-Carlo simulation creating a set of samples, return
             * energy */

            // reinitialize samples
            initializeSample();

            // save first part of progress bar
            std::string progressPosition, progressBuffer;
            if (progressDivider) {
                progressPosition = Methods::stringPos(m_rank, 3) + "Progress: "
                    "[";
            } // end if

            // make sample
            bool accepted;
            unsigned int iterations = m_maxIterations + thermalizationIdx;
            for (unsigned int i = 0; i < iterations; ++i) {
                /* sample values */
                accepted = sample();

                // gather expectation values
                accumulateSample(i, accepted);

                // print progress
                if (progressDivider) {
                    /* show progress if given */
                    if (!(static_cast<int>(fmod(i, Methods::divider(i, iterations,
                                            progressDivider))))) {
                        /* print only a few times */
                        progressBuffer = progressPosition;
                        Methods::printProgressBar(progressBuffer,
                                (float)((i==iterations-1) ? i : (i+1)) /
                                iterations, 23, "VMC");
                    } // end if
                } // end if
            } // end fori

            if (m_makeOnebodyDensities) {
                /* dump histrogram to file */
                Eigen::ArrayXd rVals = Eigen::ArrayXd::LinSpaced(m_bins,
                        m_rStep, m_rBinsMax);
                m_histogram *= wf->getNumberOfParticles() /
                    (m_histogram.sum()*rVals.pow(m_dim-1));
                writeHistogramToFile();
            } // end if

            if (m_makeRadialDensities) {
                Eigen::ArrayXd rVals = Eigen::ArrayXd::LinSpaced(m_bins,
                        m_rStep, m_rBinsMax);
                for (unsigned int i = 0; i < m_radial2DHistogram.rows(); ++i) {
                    m_radial2DHistogram.row(i) /= (m_radial2DHistogram.sum()*rVals*rVals);
                } // end fori
                writeRadialHistogramToFile();
            } // end fi

            // average expectation values
            m_accumulativeValues /= m_maxIterations;

            derivativeVariational();

            return m_accumulativeValues.energy;
        } // end function sampler

        const Eigen::VectorXd& getParameters() const {
            /* return variational parameters as a vector */
            return wf->m_parameters;
        } // end function getParameters

        const Eigen::VectorXd& getOldParameters() const {
            /* return variational parameters as a vector */
            return m_oldParameters;
        } // end function getParameters

        const Eigen::VectorXd& getNewDerivativeParameters() const {
            /* return new derivative vector */
            return m_newDerivativeParameters;
        } // end function getDerivative

        const Eigen::VectorXd& getOldDerivativeParameters() const {
            /* return new derivative vector */
            return m_oldDerivativeParameters;
        } // end function getDerivative

        const double& getAcceptance() const {
            /* return acceptance rate */
            return m_accumulativeValues.acceptance; 
        } // end function getEnergySquared

        const double& getEnergy() const {
            /* return sampled energy */
            return m_accumulativeValues.energy;
        } // end function getEnergy

        const double& getEnergySquared() const {
            /* return sampled squared energy */
            return m_accumulativeValues.energySquared; 
        } // end function getEnergySquared

        const double& getPotentialEnergy() const {
            /* return sampled potential energy */
            return m_accumulativeValues.potentialEnergy; 
        } // end function getPotentialEnergy 

        const double& getKineticEnergy() const {
            /* return sampled kinetic energy */
            return m_accumulativeValues.kineticEnergy; 
        } // end function getKineticEnergy 

        void setOneBodyDensities(std::string filename, unsigned int nBins,
                double nRmax) {
            /* switch on accumulating one-body densities */
            m_oneBodyDensityFilename = filename;
            m_bins = nBins;
            m_rBinsMax = nRmax;
            m_rStep = m_rBinsMax / m_bins;
            m_histogram = Eigen::ArrayXd::Zero(nBins);
            m_makeOnebodyDensities = true;
        } // end function setOneBodyDensities

        void setRadialDensities(std::string filename, unsigned int nBins,
                double nRmax) {
            m_radialDensityFilename = filename;
            m_bins = nBins;
            m_rBinsMax = nRmax;
            m_rStep = m_rBinsMax / m_bins;
            m_radial2DHistogram = Eigen::ArrayXXd::Zero(nBins, nBins);
            m_makeRadialDensities = true;
        } // end function setRadialDensities

        void setParameters(Eigen::VectorXd parameters) {
            /* set parameters */
            if (parameters.size() != numParameters) {
                numParameters = parameters.size();
            } // end if
//             for (unsigned int i = 0; i < parameters.size(); ++i) {
//                 /* make sure parameters are positive */
//                 if (parameters(i) <= 0) {
//                     parameters(i) = 1e-8;
//                 } // end if
//             } // end fori
            wf->setParameters(parameters);
        } // end function setParameters

        void setStep(const double &value) {
            /* set step-length */
            stepMonteCarlo = value;
            sqrtStep = sqrt(value);
        } // end function setStep

        void setMaxIterations(const unsigned int& m) {
            /* set max number of Monte Carlo samples */
            m_maxIterations = m;
        } // end function setMaxIterations

        void setNumParticles() {
            /* grab number of particles from wavefunction */
            numParticles = wf->getNumberOfParticles();
            discreteDistribution = std::uniform_int_distribution<unsigned int>(0,
                    numParticles-1);
        } // end function setNumParticles

        void setResampler(Resampler *r) {
            /* set resampler object */
            rs = r;
        } // end function setResamplingMethod

    protected:
        Accumulative m_accumulativeValues;

        double sqrtStep, stepMonteCarlo;
        bool m_hasSampleSetup;
        unsigned int numParticles, numParameters, dIdx, m_maxIterations, m_dim;
        
        std::mt19937_64 genMT64;
        std::uniform_real_distribution<double> realDistribution;

        double MetropolisRatio(double pIdx) {
            /* return ratio between new and old distributions (wavefunctions) */
            return wf->wavefunctionRatio(pIdx);
        } // end function metropolisRatio

        bool MetropolisTest(double ratio) {
            return ((ratio >= realDistribution(genMT64)) ? true : false);
        } // end function MetropolisTest

        unsigned int randomParticle() {
            return discreteDistribution(genMT64);
        } // end function randomParticle

        virtual void initializeSample() = 0;
        virtual bool sample() = 0;
};

#endif /* VMC_H */
