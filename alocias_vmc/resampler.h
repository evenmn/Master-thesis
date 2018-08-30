#ifndef RESAMPLER_H
#define RESAMPLER_H

#include <string>
#include <eigen3/Eigen/Dense>
#include <mpi.h>
#include <random>

class Resampler {
    private:
        void bootStrap(unsigned int);
        void blocking(unsigned int);
        void autoBlocking(unsigned int);

        void (Resampler::*resamplingFunction)(unsigned int);

        void distributor(const unsigned int&, unsigned int&, int*&);
        void gatherer(const Eigen::ArrayXd&, const int*, const unsigned int&,
                Eigen::ArrayXd&);
       
        std::mt19937_64 genMT64;
        std::uniform_int_distribution<int> discreteDistribution;

        int m_rank, m_procs;

        MPI_Comm comm;

    public:
        Resampler (const std::string&, MPI_Comm);
        virtual ~Resampler ();
        
        double mean, variance, std;

        Eigen::ArrayXd samples;
        
        Eigen::ArrayXd blockSamples;
        Eigen::ArrayXd blockSizes;

        void setMethod(const std::string&);
        void setSampleSize(const unsigned int&);
        void setSample(const Eigen::ArrayXd&);

        void runResampling(const unsigned int&);

        void calculateMean();
};

#endif /* RESAMPLER_H */
