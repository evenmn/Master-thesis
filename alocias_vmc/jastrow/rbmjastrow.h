#ifndef RBMJASTROW_H
#define RBMJASTROW_H

#include <vector>

#include <Eigen/Dense>

class SlaterJastrow;

class RBMJastrow {
    private:
        unsigned int m_parametersDispl, m_numJastrowParameters,
                     m_numHiddenBias, m_numVisibleBias, m_numWeights,
                     m_hiddenDispl, m_weightsDispl;
        
        Eigen::VectorXd m_jastrowGradientVector;

        Eigen::VectorXd m_expSumVector, m_oldExpSumVector;

        Eigen::MatrixXd m_weightsSumMatrix, m_oldWeightsSumMatrix;
        
        SlaterJastrow* SJ;

        void setWeightsSum();

        unsigned int wi(const unsigned int&,const unsigned int&, const unsigned
                int&);

    public:
        RBMJastrow (SlaterJastrow*);
        virtual ~RBMJastrow ();

        void resetGradient(const unsigned int&);
        void acceptGradient(const unsigned int&);
        void calculateGradient();
        void calculateGradient(const unsigned int&);

    protected:
        std::vector<double(RBMJastrow::*)(const unsigned int&)>
            jastrowVariationalDerivativeFunctionList;
        
        double jastrowWavefunctionRatio();
        double jastrowWavefunctionRatio(const unsigned int&);
        double derivativeVisibleBias(const unsigned int&);
        double derivativeHiddenBias(const unsigned int&);
        double derivativeWeight(const unsigned int&);
        double jastrowLaplacian();

        const Eigen::VectorXd& gradient(const unsigned int&);

        void initializeMatrices(const unsigned int&);
};

#endif /* RBMJASTROW_H */
