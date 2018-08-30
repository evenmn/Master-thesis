#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "vmc.h"
#include "methods.h"
#include "MTLS.h"

template<class T> class VMC;

template<class T>
class Minimizer {
    private:
        VMC<T>* vmc;

        double A, fmax, fmin, a, w, step, tol, maxStepASGD, prevValue, mu, eta,
               eps, exPol, delta, bestEnergy, bestVariance, maxValue,
               oldEnergy, temp, beta, astep, annealingFraction, threshAccept,
               tempMax;

        struct ParamsMTLS {
            /* struct of default parameters for line search in CG */
            double maxIterations;
            double mu;
            double eta;
            double delta;
            double bisectWidth;
            double bracketTol;
            double aMin0;
            double aMax0;
        } pMTLS; // end struct p1

        unsigned int m_runMax, maxIterationsLinesearch, prevIdx, sianIdx;

        bool hasSetup;

        std::string method, currMethod;

        Eigen::ArrayXd tVec;

        Eigen::VectorXd searchDirection, bestParameters;

        Eigen::MatrixXd hessianInverse, hessian;

        bool FTSIP() {
            /* Fuck-This-Sample-In-Particular. Check if sample fucked up(was
             * unstable) and fuck it (straight up) if it did, aka find one of
             * its neighbours (who hopefully isnt as shitty) and go for that
             * one. If it didnt fuck up do nothing and calmly cheer for the
             * sample and thank the fuck for it keeping it real. */
            if ((vmc->getEnergy() < 0) ||
                    (fabs(vmc->getNewDerivativeParameters().norm()) > 1000) ||
                    (vmc->getAcceptance() < threshAccept) ||
                    (sqrt((vmc->getEnergySquared() - vmc->getEnergy() *
                           vmc->getEnergy()) / vmc->m_maxIterations > 1.0))) {
                /* slightly shift and return false */
                vmc->setParameters(vmc->m_oldParameters . 
                unaryExpr([](double val) {
                        std::normal_distribution<double> nd(val, 1e-4);
                        static std::mt19937_64 rng(std::stoi(std::to_string(
                                        std::chrono::
                                        high_resolution_clock::now() .
                                        time_since_epoch() .
                                        count()).substr(10)));
                        return nd(rng);
                    }));
                vmc->sampler();
                return false;
            } else {
                return true;
            } // end if
        } // end function FTSIP

        void (Minimizer::*minimizeFunction)();

        void setParamsASGD() {
            /* set specific parameters used in ASGD method */
            A = 60.0;
            fmax = 1.0;
            fmin = -0.5;
            a = 0.9;
            w = 0.5;
            maxStepASGD = 0.1;

            tVec = Eigen::ArrayXd::Constant(vmc->numParameters, A);
        } // end function setParamsASGD

        void setParamsMTLS() {
            /* set parameters used in line search in CG and BFGS method */
            pMTLS.maxIterations = 10;
            pMTLS.mu = 0.001;
            pMTLS.eta = 0.6;
            pMTLS.delta = 4.0;
            pMTLS.bisectWidth = 0.66;
            pMTLS.bracketTol = 1e-14;
            pMTLS.aMin0 = 0.0;
            pMTLS.aMax0 = 2.0;
        } // end function setParamsMTLS

        void setParamsSABFGS() {
            /* set parameters used in SABFGS method */
            beta = 0.6;
        } // end function setParamsSABFGS
        
        void minimizeSD() {
            /* find optimal variational parameters */
            if (!hasSetup) {
                step = 0.5;

                hasSetup = true;
            } // end if

            // find new values
            Methods::steepestDescent<double>(step, vmc->wf->m_parameters,
                    static_cast<Eigen::VectorXd>(vmc->m_newDerivativeParameters
                        / vmc->m_newDerivativeParameters.norm()));
            if (FTSIP()) {
                vmc->sampler();
            } // end if
        } // end function minimizeSD

        void minimizeSABFGS() {
            /* minimize with Stochastic-Adaptive-BFGS method */
            if (!hasSetup) {
                /* allocate for minimization */
                setParamsSABFGS();

                hessianInverse = Eigen::MatrixXd::Identity(vmc->numParameters,
                        vmc->numParameters);
                hessian = Eigen::MatrixXd::Identity(vmc->numParameters,
                        vmc->numParameters);
                searchDirection = Eigen::VectorXd::Zero(vmc->numParameters);

                hasSetup = true;
            } // end if

            static double incrementVal = 1.01;

            auto setSteps = [this]() {
                /* set parameters used in update */
                delta = sqrt(Methods::innerProd(searchDirection, hessian));
                astep = Methods::innerProd(vmc->m_newDerivativeParameters,
                        hessianInverse) / (delta*delta);
                step = 1 / (1./astep + delta);
            }; // end lambda setSteps
            
            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;

            // set search direction (using derivative from previous iteration)
            searchDirection = -hessianInverse * vmc->m_oldDerivativeParameters;

            /* set parameters used in update */
            setSteps();
            
            // sample new values;
            vmc->setParameters(vmc->m_oldParameters + step*searchDirection);
            vmc->sampler();
            if (FTSIP()) {
                /* make sure sample is good */
                if (Methods::strongWolfeCondition(beta,
                            vmc->m_newDerivativeParameters,
                            vmc->m_oldDerivativeParameters, searchDirection)) {
                    /* Make Stochastic-Adaptive step if Strong Wolfe-Conditions are
                     * met */
                    searchDirection = -vmc->m_oldDerivativeParameters;
                    setSteps();
                } else {
                    /* Update Hessian and its inverse with BFGS */
                    Methods::BFGSInverse<double>(hessianInverse,
                            vmc->wf->m_parameters, vmc->m_oldParameters,
                            vmc->m_newDerivativeParameters,
                            vmc->m_oldDerivativeParameters);
                    Methods::BFGS<double>(hessian, vmc->wf->m_parameters,
                            vmc->m_oldParameters, vmc->m_newDerivativeParameters,
                            vmc->m_oldDerivativeParameters);
                } // end ifelse

                // update parameters
                vmc->setParameters(vmc->m_oldParameters + step*searchDirection);
            } else {
                /* reset hessians if sample is not good */
                hessianInverse =
                    Eigen::MatrixXd::Identity(hessianInverse.rows(),
                            hessianInverse.cols());
                hessian = Eigen::MatrixXd::Identity(hessian.rows(),
                        hessian.cols());
            } // endifelse
        } // end function minimizeSABFGS

        void minimizeBFGS() {
            /* find optimal variational parameters */
            if (!hasSetup) {
                /* allocate for minimization */
                setParamsMTLS();

                hessianInverse = Eigen::MatrixXd::Identity(vmc->numParameters,
                        vmc->numParameters);
                searchDirection = Eigen::VectorXd::Zero(vmc->numParameters);

                hasSetup = true;
            } // end if
            
            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;
           
            // set search direction
            searchDirection = - hessianInverse *
                vmc->m_oldDerivativeParameters;
            double sdNorm = searchDirection.norm();
            if (sdNorm >= 1e-10) {
                /* dont normalize if close to minimum for stability */
                searchDirection /= sdNorm;
            } // end if

            // run linesearch
            double s = MTLS::linesearchMoreThuente<>(&pMTLS, searchDirection,
                    vmc->m_oldParameters, vmc->m_accumulativeValues.energy,
                    vmc, static_cast<double(VMC<T>::*)(const Eigen::VectorXd&,
                        const unsigned int)>(&VMC<T>::sampler),
                    &VMC<T>::getNewDerivativeParameters, 0);

            // sample with new parameter and set derivates and update hessian
            // (inverse) with Broyden-Fletcher-Goldfarb-Shanno method
            vmc->setParameters(vmc->m_oldParameters + s*searchDirection);
            vmc->sampler();
            if (FTSIP()) {
                Methods::BFGSInverse<double>(hessianInverse,
                        vmc->wf->m_parameters, vmc->m_oldParameters,
                        vmc->m_newDerivativeParameters,
                        vmc->m_oldDerivativeParameters);
            } else {
                /* start anew if we FTSIPed */
                hessianInverse =
                    Eigen::MatrixXd::Identity(hessianInverse.rows(),
                            hessianInverse.cols());
            } // end ifselse
        } // end function minimizeBFGS

        void minimizeCG() {
            /* find optimal variational parameters with the Polak-Ribiere
             * method method */
            if (!hasSetup) {
                /* initialize */
                setParamsMTLS();
                searchDirection = -vmc->m_newDerivativeParameters;

                hasSetup = true;
            } // end if

            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;
           
            // perform linesearch and update parameters
            double s = MTLS::linesearchMoreThuente<>(&pMTLS, searchDirection,
                    vmc->m_oldParameters, vmc->m_accumulativeValues.energy,
                    vmc, static_cast<double(VMC<T>::*)(const Eigen::VectorXd&,
                        const unsigned int)>(&VMC<T>::sampler),
                    &VMC<T>::getNewDerivativeParameters, 0);

            // evaluate new derivatives and values and update search direction
            // according to Polak-Ribiere method (forcing b to be such that the
            // search direction is always a descent direction)
            vmc->setParameters(vmc->m_oldParameters + s*searchDirection);
            vmc->sampler();
            if (FTSIP()) {
                /* Only update is sample behaves */
                double b =
                    Methods::max((vmc->m_newDerivativeParameters.squaredNorm() -
                                vmc->m_newDerivativeParameters.transpose() *
                                vmc->m_oldDerivativeParameters) /
                            vmc->m_oldDerivativeParameters.squaredNorm(), 0.0);
                searchDirection = b * searchDirection -
                    vmc->m_newDerivativeParameters;
            } // end if
        } // end function minimizeCG

        void minimizeASGD() {
            /* find optimal variational parameters with Adaptive Stochastic
             * Gradient Descent */
            if (!hasSetup) {
                /* initialize */
                setParamsASGD();

                hasSetup = true;
            } // end if

            // array of values with damping step-function
            Eigen::ArrayXd f = fmin + (fmax - fmin) / (1 - (fmax/fmin) *
                    exp(-(vmc->m_newDerivativeParameters.array() *
                            vmc->m_oldDerivativeParameters.array() / w)));

            // set variables parameters involved in step
            tVec += f;
            for (unsigned int i = 0; i < tVec.size(); ++i) {
                if (tVec(i) < 0) {
                    tVec(i) = 0;
                } // end if
            } // end fori

            // set steps
            Eigen::ArrayXd stepVec = a / (tVec + A);

            // limit to a max step value
            for (unsigned int i = 0; i < tVec.size(); ++i) {
                if (fabs(stepVec(i)) > maxStepASGD) {
                    stepVec(i) *= maxStepASGD / fabs(stepVec(i));
                } // end if
            } // end fori

            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;

            // update parameters and derivatives
            vmc->setParameters(vmc->wf->m_parameters - (stepVec *
                        vmc->m_newDerivativeParameters.array() /
                        vmc->m_newDerivativeParameters.norm()).matrix());
            vmc->sampler();
            FTSIP();
        } // end function minimizeASGD

        void minimizeSIAN() {
            /* minimize with simulated annealing */
            if (!hasSetup) {
                /* initialize */
                vmc->m_oldDerivativeParameters =
                    vmc->m_newDerivativeParameters;
                bestParameters = vmc->wf->m_parameters;
                bestEnergy = vmc->getEnergy();
                bestVariance = (vmc->getEnergySquared() -
                        bestEnergy*bestEnergy) / vmc->m_maxIterations;
                oldEnergy = vmc->getEnergy();
                maxValue = 10.;
                tempMax = 50.;
                temp = tempMax;
                sianIdx = 1;

                hasSetup = true;
            } // end if

            // propose new state in a unit gaussian around parameter
            Eigen::ArrayXd newParameters =
                Eigen::ArrayXd::Zero(vmc->wf->m_parameters.size()) .
                unaryExpr([](double val) {
                        std::normal_distribution<double> nd(val, 1.0);
                        static std::mt19937_64 rng(std::stoi(std::to_string(
                                        std::chrono::
                                        high_resolution_clock::now() .
                                        time_since_epoch() .
                                        count()).substr(10)));
                        return nd(rng);
                    });
            double energy = vmc->sampler(newParameters);
            if (!FTSIP()) {
                /* dont bother with updating if sample is shit */
                return;
            } // end if

            // use Metropolis-Test to accept/reject state based on the gradient
            // ratio
            double energyDiff = oldEnergy - energy;
            if (vmc->MetropolisTest(exp(energyDiff/temp))) {
                /* update parameters (accept new state) */
                vmc->m_oldParameters = vmc->wf->m_parameters;
                vmc->m_oldDerivativeParameters =
                    vmc->m_newDerivativeParameters;
                double tmpVarience = (vmc->getEnergySquared() - energy*energy)
                    / vmc->m_maxIterations;
                oldEnergy = energy;
                if ((tmpVarience < bestVariance) && (energy < bestEnergy) &&
                        (vmc->getAcceptance() > 0.4)) {
                    /* make choice based on energy variance */
                    bestParameters = vmc->wf->m_parameters;
                    bestEnergy = energy;
                    bestVariance = tmpVarience;
                } // end if
            } else {
                /* reset */
                vmc->setParameters(vmc->m_oldParameters);
                vmc->m_newDerivativeParameters =
                    vmc->m_oldDerivativeParameters;
            } // end if

            temp = tempMax/sianIdx;
            sianIdx++;
        } // end function minimizeSian

        void setup() {
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters =
                Eigen::VectorXd::Constant(vmc->numParameters, 1000);
            vmc->m_newDerivativeParameters =
                Eigen::VectorXd::Zero(vmc->numParameters);
        } // end function setup

        void showProgress(const unsigned int &index) {
            /* show progress bar for each process */

            // set progress string and buffer
            const std::string progressPosition = Methods::stringPos(vmc->m_rank, 3)
                + "Progress: [";
            std::string progressBuffer;

            // value determining the frequency of progress updates
            static const double divValue = exp(fmod(4,m_runMax));
            if (!static_cast<int>(fmod(index, Methods::divider(index, m_runMax,
                                divValue)))) {
                progressBuffer = progressPosition;
                Methods::printProgressBar(progressBuffer, (float)(((index==m_runMax-1)
                                ? index : index+1)) / m_runMax, 50, currMethod);
            } // end if
        } // end function showProgress

        inline bool updater(const unsigned int& m) {
            /* wrapper for updating minimization methods in function minimize,
             * return false in case break threshold is reached and false if not
             * (return true by default) */

//             double testNorm = vmc->m_newDerivativeParameters.norm() /
//                 Methods::max(1.0, vmc->wf->m_parameters.norm());
            double testNorm = vmc->m_newDerivativeParameters.norm();
            if (testNorm <= tol) {
                /* break if values are sufficiently convergent */
                return false;
            } else {
                /* check for switches */
                if (!currMethod.compare("SIAN") && (m > static_cast<unsigned
                            int>(ceil(m_runMax*annealingFraction) + 1)) &&
                        method.compare("SIAN")) {
                    /* Stop annealing when exhausted (unless main method is
                     * SIAN) and use current best parameters */
                    vmc->setParameters(bestParameters);
                    if (!method.compare("BFGSADA") || !method.compare("BFGS"))
                    {
                        /* switch to BFGS */
                        currMethod = "BFGS";
                        minimizeFunction = &Minimizer::minimizeBFGS;
                    } else if (!method.compare("SABFGS")) {
                        /* switch to SABFGS */
                        currMethod = "SABFGS";
                        minimizeFunction = &Minimizer::minimizeSABFGS;
                    } else if (!method.compare("CGADA") ||
                            !method.compare("CG")) {
                        /* switch to CG */
                        currMethod = "CG";
                        minimizeFunction = &Minimizer::minimizeCG;
                    } else {
                        /* switch to ASGD */
                        currMethod = "ASGD";
                        minimizeFunction = &Minimizer::minimizeASGD;
                    } // end ifeif
                    hasSetup = false;
                } // end if

                if ((testNorm <= 10*tol) && (!method.compare("BFGSADA")
                            || !method.compare("CGADA"))) {
                    /* change to adaptive stochastic gradient descent method if
                     * sufficiently close to minimum */
                    currMethod = "ASGD";
                    minimizeFunction = &Minimizer::minimizeASGD;
                    hasSetup = false;
                } // end if
            } // end ifelse

            if (!method.compare("SD")) {
                /* reduce step size when close to minimum for steepest descent
                 * */
                if (testNorm <= 100*tol) {
                    step = 0.005;
                } else if (testNorm <= 10*tol) {
                    step = 0.001;
                } // end ifeif
            } // end ifeifelse

            return true;
        } // end function updater
    
    public:
        Minimizer(VMC<T> *vIn, std::string minimizationMethod, unsigned int
                maxMinimizationSamples, double threshold) {
            /* set vmc object */
            vmc = vIn;

            m_runMax = Methods::divider(vmc->m_rank, maxMinimizationSamples,
                    vmc->m_numprocs);

            tol = threshold;

            method = minimizationMethod;
            currMethod = minimizationMethod;
            prevIdx = 10;

            hasSetup = false;

            // always anneal before running BFGS, CG or ASGD
//             currMethod = "SABFGS";
//             minimizeFunction = &Minimizer::minimizeSABFGS;
            currMethod = "SIAN";
            minimizeFunction = &Minimizer::minimizeSIAN;
            annealingFraction = 0.1;

            // set function pointer to specific minimize method
            if (!minimizationMethod.compare("SD")) {
                /* Steepest Descent */
                currMethod = "SD";
                minimizeFunction = &Minimizer::minimizeSD;
            } // end if

            checkMethod(minimizationMethod);
        } // end constructor
        
        virtual ~Minimizer () {
        } // end deconstructor

        void checkMethod(const std::string& meth) {
            /* make sure method is implemented */
            if (!(!meth.compare("SD") ||
                  !meth.compare("CG") ||
                  !meth.compare("BFGS") ||
                  !meth.compare("CGADA") ||
                  !meth.compare("BFGSADA") ||
                  !meth.compare("SABFGS") ||
                  !meth.compare("ASGD") ||
                  !meth.compare("SIAN"))) {
                /* throw error if none of the above is given */
                throw std::runtime_error("Please give method, SD, BFGS, ASGD "
                        "BFGSADA, CG or CGADA");
            } // end if
        } // end function checkMethod

        void setMethod(const std::string& meth) {
            /* set method to be used */
            checkMethod(meth);
            method = meth;
            currMethod = meth;
            annealingFraction = 0.0;
            hasSetup = false;
            if (!meth.compare("SD")) {
                minimizeFunction = &Minimizer::minimizeSD;
            } else if (!meth.compare("CG") || !meth.compare("CGADA")) {
                minimizeFunction = &Minimizer::minimizeCG;
            } else if (!meth.compare("BFGS") || !meth.compare("BFGSADA")) {
                minimizeFunction = &Minimizer::minimizeBFGS;
            } else if (!meth.compare("SABFGS")) {
                minimizeFunction = &Minimizer::minimizeSABFGS;
            } else if (!meth.compare("ASGD")) {
                minimizeFunction = &Minimizer::minimizeASGD;
            } else if (!meth.compare("SIAN")) {
                minimizeFunction = &Minimizer::minimizeSIAN;
            } else {
                throw std::runtime_error("Please give method, SD, BFGS, ASGD "
                        "BFGSADA, CG or CGADA");
            } // end ifeif---else

        } // end function setMethod

        void setAnnealingFraction(double a) {
            /* set the fraction at which annealing should proceed */
            annealingFraction = a;
        } // end function setAnnealingFraction

        void minimize(bool showBar=true) {
            /* sample and find optimal parameters */
            
            // show progressbar at 0
            if (showBar) {
                showProgress(0);
            } // end if

            // initialize
            setup();
            vmc->sampler();
            prevValue = vmc->m_accumulativeValues.energy;
            std::ofstream outfile("w1.0_D2_N6_Energies.txt",
                    std::ios_base::app);
            outfile.precision(16);
            outfile.setf(std::ios::fixed);
            outfile.setf(std::ios::showpoint);
            for (unsigned int m = 0; m < m_runMax; ++m) {
                /* run minimization */

                // show progress
                if (showBar) {
                    showProgress(m);
                } // end if

                // check current state of values and update accordingly
                if (!updater(m)) {
                    break;
                } // end if
               
                // update old value
                if ((m>=prevIdx) && ((m%prevIdx)==0)) {
                    prevValue = vmc->m_accumulativeValues.energy;
                } // end if

                if (currMethod.compare("SIAN")) {
                    outfile << vmc->m_accumulativeValues.energy << " " <<
                        (vmc->m_accumulativeValues.energySquared -
                         vmc->m_accumulativeValues.energy *
                         vmc->m_accumulativeValues.energy) /
                        vmc->m_maxIterations << "\n";
                    outfile.flush();
                } // end if

                // sample and minimize with specified method
                (this->*minimizeFunction)();

                std::cout << std::setprecision(14) <<
                    Methods::stringPos(vmc->m_rank+3, vmc->m_numprocs) <<
                    vmc->getAcceptance() << "       " <<
                    vmc->m_newDerivativeParameters.norm() << "      \n\n" <<
//                     vmc->m_newDerivativeParameters.transpose() << "     \n\n" <<
//                     vmc->wf->m_parameters << "   \n\n" <<
                    "    " << vmc->m_accumulativeValues.energy << "   " <<
                    (vmc->m_accumulativeValues.energySquared -
                     vmc->m_accumulativeValues.energy *
                     vmc->m_accumulativeValues.energy) / vmc->m_maxIterations
                    << "\n\n";
            } // end form

            // show 100% in end
            if (showBar) {
                showProgress(m_runMax-1);
            } // end if
           
            if (outfile.is_open()) {
                outfile.close();
            } // end if
        } // end function minimize

        void setFunctionThresh(double thresh) {
            /* set threshold for acceptance */
            threshAccept = thresh;
        } // end function setFunctionThresh

        const Eigen::VectorXd& getDerivative() const {
            /* return new derivative vector */
            return vmc->getNewDerivativeParameters();
        } // end function getDerivative

        const Eigen::MatrixXd& getHessianInverse() const {
            /* return new derivative vector */
            return hessianInverse;
        } // end function getDerivative
};

#endif /* MINIMIZER_H */
