#ifndef PADENQS_H
#define PADENQS_H

#include <vector>

#include <Eigen/Dense>

#include "padejastrow.h"
#include "rbmjastrow.h"

class SlaterJastrow;

class PadeNQS : public PadeJastrow, public RBMJastrow {
    private:
        Eigen::VectorXd m_jastrowGradientVector;
        
        SlaterJastrow* SJ;

    protected:
        std::vector<double(PadeNQS::*)(const unsigned int&)>
            jastrowVariationalDerivativeFunctionList;
        
        const Eigen::VectorXd& gradient(const unsigned int&);

        double jastrowWavefunctionRatio();
        double jastrowWavefunctionRatio(const unsigned int&);
        double jastrowLaplacian();

        void initializeMatrices(const unsigned int&);

    public:
        PadeNQS (SlaterJastrow*);
        virtual ~PadeNQS ();
        
        void resetGradient(const unsigned int&);
        void acceptGradient(const unsigned int&);

        void calculateGradient();
        void calculateGradient(const unsigned int&);
};

#endif /* PADENQS_H */
