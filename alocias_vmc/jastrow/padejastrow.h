#ifndef PADEJASTROW_H
#define PADEJASTROW_H

#include <vector>

#include <Eigen/Dense>

class SlaterJastrow;

class PadeJastrow {
    private:
        unsigned int m_parametersDispl, m_numJastrowParameters;
        double parallelSpinFactor, antiParallelSpinFactor, dimMinusOne;

        Eigen::VectorXd m_jastrowGradientVector;
        
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
            m_jastrowGradient3DMatrix, m_oldJastrowGradient3DMatrix;
        
        double jastrowWavefunctionRatioExpression(const unsigned int&, const
                unsigned int&);
        
        const Eigen::VectorXd& jastrowGradientExpression(const unsigned int&,
                const unsigned int&);

        SlaterJastrow* SJ;

    public:
        PadeJastrow (SlaterJastrow*);
        virtual ~PadeJastrow ();

        void resetGradient(const unsigned int&);
        void acceptGradient(const unsigned int&);

        void calculateGradient();
        void calculateGradient(const unsigned int&);

    protected:
        std::vector<double(PadeJastrow::*)(const unsigned int&)>
            jastrowVariationalDerivativeFunctionList;
        
        double spinFactor(const unsigned int&, const unsigned int&);

        const Eigen::VectorXd& gradient(const unsigned int&);

        double jastrowWavefunctionRatio();
        double jastrowWavefunctionRatio(const unsigned int&);
        double jastrowVariationalDerivativeExpression(const unsigned int&);
        double jastrowLaplacian();

        void initializeMatrices(const unsigned int&);
};

#endif /* PADEJASTROW_H */
