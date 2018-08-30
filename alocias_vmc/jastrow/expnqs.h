#ifndef EXPNQS_H
#define EXPNQS_H

#include <Eigen/Dense>

#include "padejastrow.h"
#include "rbmjastrow.h"

class SlaterJastrow;

class ExpNQS : private PadeJastrow, private RBMJastrow {
    private:
        Eigen::VectorXd m_jastrowGradientVector;
        
        SlaterJastrow* SJ;
        
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
            m_jastrowGradient3DMatrix, m_oldJastrowGradient3DMatrix;
        
        double jastrowWavefunctionRatioExpression(const unsigned int&, const
                unsigned int&);
        
        const Eigen::VectorXd& jastrowGradientExpression(const unsigned int&,
                const unsigned int&);
    
    protected:
        std::vector<double(ExpNQS::*)(const unsigned int&)>
            jastrowVariationalDerivativeFunctionList;
        
        const Eigen::VectorXd& gradient(const unsigned int&);

        double jastrowWavefunctionRatio();
        double jastrowWavefunctionRatio(const unsigned int&);
        double jastrowLaplacian();

        void initializeMatrices(const unsigned int&);

    public:
        ExpNQS (SlaterJastrow*);
        virtual ~ExpNQS ();
        
        void resetGradient(const unsigned int&);
        void acceptGradient(const unsigned int&);

        void calculateGradient();
        void calculateGradient(const unsigned int&);
};

#endif /* EXPNQS_H */
