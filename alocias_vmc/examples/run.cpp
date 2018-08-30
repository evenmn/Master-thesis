#include <vector>
#include "run.h"

#ifdef HARMONICOSCILLATOR
bool checkInputFile(YAML::Node& inputs) {
    /* check that non-optional parameters are given in input file */
    if (!inputs["omega"] || 
        !inputs["numparticles"] || 
        !inputs["stepmc"] || 
        !inputs["dim"] || 
        !inputs["maxitermc"]) {
        return false;
    } else {
        return true;
    } //  end ifelse
} // end function checkInputFile
#endif

#if defined HARTREEFOCK || VARIATIONALHARTREEFOCK || HARTREEFOCKDOUBLEWELL || defined VARIATIONALHARTREEFOCKDOUBLEWELL
bool checkInputFile(YAML::Node& inputs) {
    if (!inputs["omega"] || 
        !inputs["numparticles"] || 
        !inputs["stepmc"] ||
        !inputs["dim"] || 
        !inputs["maxitermc"] || 
        !inputs["numbasis"] ||
        !inputs["coeffs"]) {
        return false;
    } else {
        return true;
    } //  end ifelse
} // end function
#endif

YAML::Node argsParser(const std::string& inputFile) {
    /* parse arguments and set non-optional parameters. return a YAML node with
     * the parsed arguments */

    // create YAML node map
    YAML::Node inputs = YAML::LoadFile(inputFile);
    if (!checkInputFile(inputs)) {
        /* check that non-optional parameters are given in input file */
        std::cout << "Input file incorrectly setup" << std::endl;
        exit(1);
    } // end if
   
    if (!inputs["importance"]) {
        /* switch on importance sampling by default */
        inputs["importance"] = true;
    } // end if
    if (!inputs["jastrow"]) {
        /* switch on Jastrow factor by default */
        inputs["jastrow"] = true;
    } // end if
    if (!inputs["progress"]) {
        /* show progress bar by default */
        inputs["progress"] = true;
    } // end if

    return inputs;
} // end function argsParser

void resample(YAML::Node& inputs, Resampler*& rs, double& resampleTime, int
        myRank) {
    /* resample in parallell */
    double myResampleTime;
    Eigen::VectorXd finalSamples =
        Eigen::VectorXd::Zero(inputs["maxitermc"].as<unsigned int>());
    if (myRank == 0) {
        finalSamples = rs->samples;
    } // end if

    #ifdef MPI_ON
    MPI_Bcast(finalSamples.data(), inputs["maxitermc"].as<unsigned int>(),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif
    rs->setSample(finalSamples);
    auto start = std::chrono::high_resolution_clock::now();
    rs->runResampling(inputs["resampling"][1].as<unsigned int>());
    auto end = std::chrono::high_resolution_clock::now();
    myResampleTime = (end - start).count();
    #ifdef MPI_ON
    MPI_Allreduce(&(myResampleTime), &(resampleTime), 1, MPI_DOUBLE,
            MPI_SUM, MPI_COMM_WORLD);
    #endif
} // end function resample

void setParameters(YAML::Node& inputs, Eigen::VectorXd& parameters, unsigned
        int K) {
    /* set initial parameters */
    if (inputs["parameters"]) {
        /* set given parameters */
        parameters = Eigen::VectorXd::Zero(inputs["parameters"].size());
        for (unsigned int i = 0; i < parameters.size(); ++i) {
            parameters(i) = inputs["parameters"][i].as<double>();
        } // end fori
    } else {
        /* default to 2 parameters */
        parameters = Eigen::VectorXd::Zero(K);
    } // end if
} // end function setParameters

double run(const std::string& inputFilename) {
    /* run sampling */

    // set inputs
    YAML::Node inputs = argsParser(inputFilename);
    double E;
    if (inputs["jastrow"].as<bool>()) {
        E = run<SlaterJastrow>(inputs);
    } else {
        E = run<Slater>(inputs);
    } // end ifselse

    return E;
} // end function run

double run(YAML::Node& inputs) {
    /* run sampling */
    double E;
    if (inputs["jastrow"].as<bool>()) {
        E = run<SlaterJastrow>(inputs);
    } else {
        E = run<Slater>(inputs);
    } // end ifselse
    
    return E;
} // end function run
