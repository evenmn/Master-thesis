#ifndef RUN_H
#define RUN_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <type_traits>
#include <cassert>

#include <eigen3/Eigen/Dense>
#include <mpi.h>
#include <yaml-cpp/yaml.h>

#include "../minimizer.h"
#include "../importanceSampling.h"
#include "../bruteforce.h"
#include "../slaterjastrow.h"
#include "../slater.h"

template<class T, class W> constexpr double threshold() {
    /* return threshold for acceptance */
    if constexpr (std::is_same_v<T, ImportanceSampling<W>>) {
        /* case importance sampling */
        return 0.989;
    } else if (std::is_same_v<T, Bruteforce<W>>) {
        /* case bruteforce */
        return 0.53;
    } else {
        /* fail if type is not Bruteforce or ImportanceSampling */
        static_assert(std::disjunction_v<std::is_same<T, Bruteforce<W>>,
                std::is_same<T, ImportanceSampling<W>>>);
    } // end if constexpr
} // end function threshold 

void resample(YAML::Node&, Resampler*&, double&, int);
void setParameters(YAML::Node&, Eigen::VectorXd&, unsigned int=1);
void runBruteForce(YAML::Node&);
void runImportanceSampling(YAML::Node&);
YAML::Node argsParser(const std::string& inputFile);
double run(const std::string&);
double run(YAML::Node& inputs);
bool checkInputFile(YAML::Node&);

#if defined HARTREEFOCKDOUBLEWELL || defined VARIATIONALHARTREEFOCKDOUBLEWELL
template<class T> T* setWavefunction(YAML::Node& inputs, const Eigen::MatrixXd&
        parameters) {
    unsigned int numParticles = inputs["numparticles"].as<unsigned int>();
    unsigned int numBasis = inputs["numbasis"].as<unsigned int>();
    YAML::Node coeffNode = inputs["coeffs"];
    Eigen::MatrixXd coeffs(numBasis, numParticles);
    for (unsigned int i = 0; i < coeffNode.size(); ++i) {
        for (unsigned int j = 0; j < coeffNode[i].size(); ++j) {
            /* stack coefficients to be the same for spin up and down(spacial
             * part is same for both spin states) */
            coeffs(j,i) = coeffNode[i][j].as<double>();
            coeffs(j,i+numParticles/2) = coeffNode[i][j].as<double>();
        } // end forj
    } // end fori

    unsigned int dim = inputs["dim"].as<unsigned int>();
    T* wf = new T(dim, numParticles, parameters);
    wf->initializeParameters(inputs["omega"].as<double>(), coeffs);
    wf->initializeMatrices(); // FIXME: Change routine...

    //FIXME: error in case Slater ir not full

    return wf;
} // end function setWavefunction
#endif

#if defined HARTREEFOCK || defined VARIATIONALHARTREEFOCK
template<class T> T* setWavefunction(YAML::Node& inputs, const Eigen::MatrixXd&
        parameters) {
    /* set wavefunction object with coefficients from file */
    unsigned int numParticles = inputs["numparticles"].as<unsigned int>();
    unsigned int numBasis = inputs["numbasis"].as<unsigned int>();
    YAML::Node coeffNode = inputs["coeffs"];
    Eigen::MatrixXd coeffs(numBasis, numParticles);
    for (unsigned int i = 0; i < coeffNode.size(); ++i) {
        for (unsigned int j = 0; j < coeffNode[i].size(); ++j) {
            /* stack coefficients to be the same for spin up and down(spacial
             * part is same for both spin states) */
            coeffs(j,i) = coeffNode[i][j].as<double>();
            coeffs(j,i+numParticles/2) = coeffNode[i][j].as<double>();
        } // end forj
    } // end fori

    unsigned int dim = inputs["dim"].as<unsigned int>();
    T* wf = new T(dim, numParticles, parameters);
    wf->initializeParameters(inputs["omega"].as<double>(), numBasis, coeffs);
    wf->initializeMatrices(); // FIXME: Change routine...

    std::string message = wf->setupDone();
    if (message.compare("")) {
        std::cout << message << std::endl;
        delete wf;
        #ifdef MPI_ON
        MPI_Finalize();
        #endif
        exit(0);
    } // end if

    return wf;
} // end function setWavefunction
#endif

#ifdef HARMONICOSCILLATOR	
template<class T> T* setWavefunction(YAML::Node& inputs, const Eigen::VectorXd&
        parameters) {
    /* set wavefunction object (assume setParameters is already called) */
    unsigned int dim = inputs["dim"].as<unsigned int>();
    unsigned int numParticles = inputs["numparticles"].as<unsigned int>();
    T* wf = new T(dim, numParticles, parameters);
    wf->initializeParameters(inputs["omega"].as<double>());
    wf->initializeMatrices(); // FIXME: Change routine...

    std::string message = wf->setupDone();
    if (message.compare("")) {
        std::cout << message << std::endl;
        delete wf;
        #ifdef MPI_ON
        MPI_Finalize();
        #endif
        exit(0);
    } // end if
    
    return wf;
} // end function setWavefunction
#endif

template<typename Sampler, typename T> void findOptimalParameters(YAML::Node&
        inputs, Sampler*& vmc, Eigen::VectorXd& initialParameters, int myRank,
        int numProcs) {
    /* run minimization with BruteForce or ImportanceSampling as Sampler */
    /* find optimal parameters */

    // set minimizer object
    Minimizer<T> *minimizer = NULL;
    if (inputs["minimization"]) {
        /* set with specified method */
        vmc->setMaxIterations(inputs["minimization"][3].as<unsigned int>());
        minimizer = new Minimizer<T>(vmc,
                inputs["minimization"][0].as<std::string>().c_str(),
                inputs["minimization"][1].as<unsigned int>(),
                inputs["minimization"][2].as<double>());

        minimizer->setFunctionThresh(threshold<Sampler, T>());
    } // end if

    // find initial parameters
    unsigned int numVariationalParameters = initialParameters.size();
    using EigenRowArrayXXd = Eigen::Array<double, Eigen::Dynamic,
          Eigen::Dynamic, Eigen::RowMajor>;
    EigenRowArrayXXd parametersBuffer = EigenRowArrayXXd::Zero(numProcs,
            numVariationalParameters);
    if (!inputs["parameters"]) {
        /* find random initial parameters if not given */
        if (myRank == 0) {
            /* let root distribute initial parameters */
            if (inputs["numhiddenbias"]) {
                /* Use Gaussian distributed for RBM */
                parametersBuffer = parametersBuffer.unaryExpr([](double dummy) {
                        static std::mt19937_64 rng(std::stoi(std::to_string(
                                        std::chrono::high_resolution_clock::now().
                                        time_since_epoch().count()).substr(10)));
                        std::normal_distribution<double> nd(0.0,0.5);
                        return nd(rng);
                    });
            } else { 
                /* default to uniform random distributed */
                double lower = 0.01;
                double higher = 2.0;
                parametersBuffer = lower + (Eigen::ArrayXXd::Random(numProcs,
                            numVariationalParameters)*0.5 +
                        0.5)*(higher-lower);
            } // end if
        } // end if
        #ifdef MPI_ON
        MPI_Scatter(parametersBuffer.data(), numVariationalParameters,
                MPI_DOUBLE, initialParameters.data(), numVariationalParameters,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
        #endif

        if (inputs["numparameters"].as<unsigned int>() > 0) {
            initialParameters.head(inputs["numparameters"].as<unsigned int>())
                = initialParameters.head(inputs["numparameters"].as<unsigned
                        int>()).unaryExpr([](double dummy) {
                        static std::mt19937_64 rng(std::stoi(std::to_string(
                                        std::chrono::high_resolution_clock ::
                                        now().time_since_epoch().count()) .
                                    substr(10)));
                        static std::uniform_real_distribution<double>
                        nd(0.9,1.7);
                        return nd(rng);
                    });
        } // end if
    } else {
        /* Dont do annealing if parameters are given */
        if (minimizer != NULL) {
            minimizer->setMethod(inputs["minimization"][0].as<std::string>());
        } // end if
    } // end if

    if ((numProcs > 1) && (minimizer != NULL)) {
        /* perform annealing fully in parallel */
        minimizer->setAnnealingFraction(1.0);
    } // end if

    vmc->setParameters(initialParameters);

    // find optimal parameters
    if (minimizer != NULL) {
        minimizer->minimize(inputs["progress"].as<bool>());
    } else {
        return;
    } // end if

    // let root gather parameters, average them and redistribute
    Eigen::ArrayXd energies = Eigen::ArrayXd::Zero(numProcs);
    Eigen::ArrayXd acceptances = Eigen::ArrayXd::Zero(numProcs);
    #ifdef MPI_ON
    MPI_Gather(&vmc->getEnergy(), 1, MPI_DOUBLE, energies.data(), 1,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&vmc->getAcceptance(), 1, MPI_DOUBLE, acceptances.data(), 1,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    parametersBuffer.setZero();
    MPI_Gather(vmc->getParameters().data(), numVariationalParameters,
            MPI_DOUBLE, parametersBuffer.data(), numVariationalParameters,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif
    if (myRank == 0 && numProcs > 1) {
        /* choose parameters which give lowest(stable) energy */
        for (unsigned int p = 0; p < numProcs; ++p) {
            if (acceptances(p) < threshold<Sampler, T>()) {
                energies(p) = 1e16;
            } // end if
        } // end forp
        vmc->setParameters(parametersBuffer.row(Methods::argMin(energies)));

        if (((numProcs > 1) && (minimizer != NULL)) ||
                inputs["minimization"][0].as<std::string>().compare("SIAN")) {
            minimizer->setMethod(inputs["minimization"][0].as<std::string>());
            minimizer->minimize(inputs["progress"].as<bool>());
        } // end if
    } // end if

    // clear memory
    delete minimizer;
} // end function findOptimalParameters

template<typename Sampler> void setResampler(YAML::Node& inputs, Sampler*& vmc,
        Resampler*& rs) {
    /* set resampler object in vmc */
    if (inputs["resampling"]) {
        rs = new Resampler(inputs["resampling"][0].as<std::string>().c_str(),
                MPI_COMM_WORLD);
        rs->setSampleSize(inputs["maxitermc"].as<unsigned int>());
        vmc->setResampler(rs);
    } // end if
} // end function setResampler

template<typename Sampler> double runFinalSampling(YAML::Node& inputs, Sampler*&
        vmc, std::chrono::duration<double>& sampleTime) {
    /* run sampling */
    auto start = std::chrono::high_resolution_clock::now();
    vmc->setStep(inputs["stepmc"].as<double>());
    double E;
    if (inputs["progress"].as<bool>()) {
        E = vmc->sampler(exp(fmod(4.5,inputs["maxitermc"].as<unsigned
                        int>())));
    } else {
        E = vmc->sampler();
    } // end if
    auto end = std::chrono::high_resolution_clock::now();
    sampleTime = end - start;

    return E;
} // end function runFinalSampling

template<typename Sampler> void finalize(YAML::Node& inputs, Sampler*& vmc,
        Resampler*& rs, std::chrono::duration<double>& sampleTime, double
        resampleTime, int numProcs) {
    /* print final values to screen or write to file if filename is given(in
     * inputs) */

    if (inputs["test"]) {
        /* dont output if testing */
        return;
    } // end if
        
    // grab values
    Eigen::VectorXd finalValues =  Eigen::VectorXd::Zero(3);
    finalValues << vmc->getAcceptance(), vmc->getEnergy(),
                vmc->getEnergySquared();

    if (inputs["output"]) {
        /* write to file if given */
        std::ofstream outfile;
        double var = (rs ?  rs->variance : (finalValues(2) -
                    finalValues(1)*finalValues(1)) /
                (inputs["maxitermc"].as<unsigned int>() - 1));
        double std = (rs ? rs->std : sqrt(var));
        outfile.open(inputs["output"].as<std::string>().c_str(),
                std::ios_base::app);
        outfile << std::setprecision(14) << "Const: " <<
            inputs["omega"].as<double>() << "\nEnergy: " << finalValues(1) <<
            "\nAcceptance: " << finalValues(0) << "\nvar: " << var << "\nstd: "
            << std << "\nParameters: " << vmc->getParameters().transpose() <<
            "\nDerivative: " <<
            vmc->getNewDerivativeParameters().transpose() << "\n" <<
            std::endl;
    } else {// end if
        /* print values */
        int inputsSize = 0;
        for (YAML::const_iterator it=inputs.begin(); it != inputs.end();
                ++it) {
            /* find number of arguments given */
            inputsSize++;
        } // end forit
        double var = (rs ?  rs->variance : (finalValues(2) -
                    finalValues(1)*finalValues(1)) /
                (inputs["maxitermc"].as<unsigned int>() - 1));
        double std = (rs ? rs->std : sqrt(var));
        inputs.remove("coeffs");
        std::cout << "\033["+std::to_string(numProcs+3)+"H" << "numprocs: "
            "" + std::to_string(numProcs) + "\n" << inputs << std::endl;
        std::cout << "\033["+std::to_string(numProcs+3)+"H" << std::endl;
        std::cout << "\033["+std::to_string(numProcs+4+inputsSize)+"H" << std::endl;
        std::cout << std::setprecision(14) << "Energy: " << finalValues(1) <<
            "\nAcceptance: " << finalValues(0) << "\nvar: " << var << "\nstd: "
            << std << "\nParameters: " << vmc->getParameters().transpose() <<
            "\nDerivative: " <<
            vmc->getNewDerivativeParameters().transpose() << "\n" <<
            "\nSampling " "duration: " << sampleTime.count() << " s" <<
            "\nResampling " "duration: " << resampleTime/numProcs << " s" <<
            std::endl;
    } // ifelse
} // end function printResults

template<class Sampler, class Wavefunction> double runSpecified(YAML::Node&
        inputs) {
    /* run with sampler Sampler */

    int myRank = 0;
    int numProcs = 1;
    #ifdef MPI_ON
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    #endif

    // set parameters if given
    Eigen::VectorXd initialParameters;
    #if defined (RBMJASTROW) || (PADENQS) || (EXPNQS)
        unsigned int numParticles = inputs["numparticles"].as<unsigned int>();
        unsigned int numHidden = inputs["numhiddenbias"].as<unsigned int>();
        unsigned int dim = inputs["dim"].as<unsigned int>();
        setParameters(inputs, initialParameters,
                inputs["numparameters"].as<unsigned int>() + numParticles*dim +
                numHidden + numParticles*dim*numHidden);
    #else
        setParameters(inputs, initialParameters,
                inputs["numparameters"].as<unsigned int>());
    #endif

    // make wavefunction and vmc object (check that Slater is full as well)
    Wavefunction* wf = setWavefunction<Wavefunction>(inputs,
            initialParameters);
    Sampler* vmc = new Sampler(wf, inputs["stepmc"].as<double>(),
            initialParameters, inputs["maxitermc"].as<unsigned int>(), myRank,
            numProcs);

    // minimize if parameters are not give
    findOptimalParameters<Sampler, Wavefunction>(inputs, vmc,
            initialParameters, myRank, numProcs);

    // initialize resampler
    Resampler *rs = NULL;
    setResampler<Sampler>(inputs, vmc, rs);

    // run final sampling on root
    std::chrono::duration<double> sampleTime;
    double E;
    if (myRank == 0) {
        if(inputs["onebody"]) {
            vmc->setOneBodyDensities(inputs["onebody"][0].as<std::string>(),
                    inputs["onebody"][1].as<unsigned int>(),
                    inputs["onebody"][2].as<double>());
        } // end if
        if(inputs["radial"]) {
            vmc->setRadialDensities(inputs["radial"][0].as<std::string>(),
                    inputs["radial"][1].as<unsigned int>(),
                    inputs["radial"][2].as<double>());
        } // end if
        vmc->setMaxIterations(inputs["maxitermc"].as<unsigned int>());
        E = runFinalSampling(inputs, vmc, sampleTime);
    } // end if

    // let root distribute whole sample to slaves if resampling is enabled
    double resampleTime = 0;
    if (rs != NULL) {
        resample(inputs, rs, resampleTime, myRank);
    } // end if


    // print final values to screen
    if (myRank == 0 && !inputs["test"]) {
        /* print only once */
        finalize<Sampler>(inputs, vmc, rs, sampleTime, resampleTime, numProcs);
    } // end if
    
    delete vmc;
    delete wf;

    // free objects
    if (rs != NULL) {
        delete rs;
    } // end if

    return E;
} // end function runSpecified

template<class Wavefunction> double run(YAML::Node& inputs) {
    /* run with wavefunction T */
    double E;
    if (inputs["importance"].as<bool>()) {
        E = runSpecified<ImportanceSampling<Wavefunction>,
          Wavefunction>(inputs);
    } else {
        E = runSpecified<Bruteforce<Wavefunction>, Wavefunction>(inputs);
    } // end if

    return E;
} // end function run

#endif /* RUN_H */
