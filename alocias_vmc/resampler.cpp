#include "resampler.h"

#include <cstring>
#include <chrono>
#include <cmath>
#include <iostream>
#include <array>

#include "methods.h"

Resampler::Resampler(const std::string& method, MPI_Comm c) {
    setMethod(method);
    
    std::mt19937_64 genMT64(std::stoi(std::to_string(
                    std::chrono::high_resolution_clock::now().
                    time_since_epoch().count()).substr(10)));

    comm = c;

    MPI_Comm_rank(comm, &m_rank);
    MPI_Comm_size(comm, &m_procs);
} // end constructor

Resampler::~Resampler() {
} // end deconstructor

void Resampler::setMethod(const std::string& method) {
    /* set method for resampling */
    if (!method.compare("bootstrap")) {
        resamplingFunction = &Resampler::bootStrap;
    } else if (!method.compare("blocking")) {
        resamplingFunction = &Resampler::blocking;
    } else if (!method.compare("autoblocking")) {
        resamplingFunction = &Resampler::autoBlocking;
    } else {
        std::cerr << "Specify resampling method" << std::endl;
        exit(1);
    } // end ifeifeif
} // end function setMethod

void Resampler::setSampleSize(const unsigned int& size) {
    /* set size of sample array */
    samples = Eigen::ArrayXd::Zero(size);
    discreteDistribution = std::uniform_int_distribution<int>(0, size-1);
} // end function setSampleSize

void Resampler::setSample(const Eigen::ArrayXd& inputSample) {
    /* set sample and sizes */
    samples = inputSample;
    discreteDistribution = std::uniform_int_distribution<int>(0,
            samples.size()-1);
} // end function setSample

void Resampler::distributor(const unsigned int &totalSize, unsigned int
        &buffer, int*& totalSizes) {
    /* let root divide totalSize into integer equal chunks and distribute each
     * size to slaves. */

    // allocate space for array if empty
    if (!totalSizes) {
        totalSizes = new int[m_procs];
    } // end if

    // allocate arrays for sizes to be distributed (only for root)
    if (m_rank == 0) {
        for (int p = 0; p < m_procs; ++p) {
            totalSizes[p] = Methods::divider(p, totalSize, m_procs);
        } // end forp
    } // end if

    // send respetive size to each slave (save in buffer)
    MPI_Scatter(totalSizes, 1, MPI_INT, &buffer, 1, MPI_INT, 0, comm);
} // end function distributor

void Resampler::gatherer(const Eigen::ArrayXd &sendBuf, const int* sizes, const
        unsigned int &mySize, Eigen::ArrayXd &buffer) {
    /* let root gather values from sendBuf into buffer */

    // temporary arrays used by gather (only relevant for root, is NULL for
    // slaves)
    int *recvc = NULL;
    int *displ = NULL;

    // let root allocate and set above arrays
    if (m_rank == 0) {
        recvc = new int[m_procs];
        displ = new int[m_procs];
        Methods::setDisplacement(recvc, displ, sizes, m_procs);
    } // end if

    // root gathers results
    MPI_Gatherv(sendBuf.data(), mySize, MPI_DOUBLE, buffer.data(), recvc,
            displ, MPI_DOUBLE, 0, comm);

    // free temporary arrays used by gather
    if (m_rank == 0) {
        delete[] recvc;
        delete[] displ;
    } // end if
} // end function gatherer

void Resampler::bootStrap(unsigned int bootSize) {
    /* resample with bootstrap method and set variance */

    // bootSizes for processes and array of sizes (only used by root, NULL for
    // slaves)
    unsigned int myBootSize;
    int* bootSizes = new int[m_procs];
    
    // root distributes chunks to slaves
    distributor(bootSize, myBootSize, bootSizes);

    // each process works on its respective chunk
    Eigen::ArrayXd bootSample = Eigen::ArrayXd::Zero(myBootSize);
    for (unsigned int i = 0; i < myBootSize; ++i) {
        double tmp = 0;
        for (unsigned int j = 0; j < samples.size(); ++j) {
            tmp += samples(discreteDistribution(genMT64));
        } // end forj
        bootSample(i) = tmp;
    } // end fori
   
    // let root gather results into total
    Eigen::ArrayXd bootSampleTotal;
    if (m_rank == 0) {
        bootSampleTotal = Eigen::ArrayXd::Zero(bootSize);
    } // end if
    gatherer(bootSample, bootSizes, myBootSize, bootSampleTotal);

    // free sizes array and calculate final mean and variance
    if (m_rank == 0) {
        delete[] bootSizes;
        bootSampleTotal /= samples.size();
        mean = bootSampleTotal.mean();
        variance = Methods::var<double>(bootSampleTotal, mean);
        std = sqrt(variance);
    } // end if
} // end function bootStrap

void Resampler::autoBlocking(unsigned int i) {
    /* resample with blocking method and set variance and covariance */
    int d = log2(samples.size());
    Eigen::ArrayXd s = Eigen::ArrayXd::Zero(d);
    Eigen::ArrayXd gamma = Eigen::ArrayXd::Zero(d);
    double mu = samples.mean();
    Eigen::ArrayXd newx = samples.array();
    unsigned int newSize = newx.size();
    Eigen::ArrayXd ev = Eigen::ArrayXd::Zero(newSize);
    Eigen::ArrayXd od = Eigen::ArrayXd::Zero(newSize);
    for (unsigned int i = 0; i < d; ++i) {
        int nm1 = newSize - 1;
        gamma(i) = 1./newSize * ((newx.segment(0,nm1) - mu) *
                (newx.segment(1,nm1) - mu)).sum();

        s(i) = Methods::var<double>(newx);

        for (unsigned int k = 0, l=0; k < newSize; k+=2, ++l) {
            od(l) = newx(k+1);
            ev(l) = newx(k);
        } // end fork
        newSize /= 2;
        newx.head(newSize) = 0.5 * (ev.head(newSize) + od.head(newSize));
    } // end fori

    Eigen::ArrayXd M = Eigen::ArrayXd::Zero(d);
    for (unsigned int k = 0; k < d; ++k) {
        M(k) = pow(2,k+1);
    } // end fork
    M *= ((gamma/s).pow(2) * M.reverse()).reverse();

    // cumulative sum
    for (unsigned int k = 0; k < d; ++k) {
        M(k) = M.head(k).sum();
    } // end fork
    M = M.reverse();

    static constexpr std::array<double, 30> q{
        6.634897, 9.210340, 11.344867, 13.276704, 15.086272, 16.811894,
        18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967,
        27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306,
        36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820,
        44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181
    };

    unsigned int k = 0;
    for (k = 0; k < d; ++k) {
        if (M(k) < q[k]) {
            break;
        } // end if
    } // end fork

    if (k >= (d-1)) {
        // FIXME: MAKE ERROR
    } // end if

    // set values
    variance = s(k)/pow(2,d-k);
    std = sqrt(variance);
} // end function autoBlocking

void Resampler::blocking(unsigned int blocks) {
    /* resample with blocking method and set variance and covariance */
    if (m_rank == 0) {
        const int minBlockSize = (int)(blocks / 2);
        const int maxBlockSize = (int)(samples.size() / minBlockSize);
        const int blockStep = (maxBlockSize - minBlockSize + 1) / blocks;
        blockSamples = Eigen::ArrayXd::Zero(blocks);
        blockSizes = Eigen::ArrayXd::Zero(blocks);
        for (unsigned int i = 0; i < blocks; ++i) {
            blockSizes(i) = maxBlockSize - i*blockStep;
            int tmpBlockNumber = (int)ceil(samples.size() /
                    ((double)blockSizes(i)));
            Eigen::ArrayXd tmpMean = Eigen::ArrayXd::Zero(tmpBlockNumber);
            for (int j = 0; j < tmpBlockNumber; ++j) {
                tmpMean(i) = samples.segment(j*blockSizes(i),
                        blockSizes(i)).mean();
            } // end forj
            blockSamples(i) = Methods::var<double>(tmpMean) / tmpBlockNumber;
        } // end fori
    } // end if
} // end function blocking 

void Resampler::runResampling(const unsigned int& s) {
    /* run resampling with specified method */
    (this->*resamplingFunction)(s);
} // end function runResampling

void Resampler::calculateMean() {
    mean = samples.mean();
} // end function
