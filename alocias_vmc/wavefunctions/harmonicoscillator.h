#ifdef HARMONICOSCILLATOR

#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H

#include "../basis/harmonicoscillatorbasis.h"

#include <Eigen/Dense>
#include <vector>

class Slater;

class Quantumdot : public QuantumdotBasis {
    private:
        double alpha, omega, omegaSq, aw, sqaw;
        
        bool m_interaction;
        
        Eigen::MatrixXd m_SnewPositions, m_SoldPositions;

        Eigen::VectorXd m_laplacianSumVec;
        
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
            m_hermite3DMatrix, m_oldHermite3DMatrix;

        void setHermite3DMatrix(const unsigned int&);

        Slater* slater;

    public:
        Quantumdot(Slater*);
        virtual ~Quantumdot();

        bool isFull;

        std::string setupDone();

        void setParameters(const Eigen::VectorXd&);
        void initializeParameters(const double&);
        void setInteraction(bool);

        double potentialEnergy();

        double calculateWavefunction(const unsigned int&, const unsigned int&);
        
        double gradientExpression(const unsigned int&, const int&, const unsigned
                int&);
        const Eigen::VectorXd& laplacianExpression(const unsigned int&, const
                unsigned int&);
        double variationalDerivativeExpression(const unsigned int&, const
                unsigned int&, const unsigned int&);

        void reSetAll();
        void initializeMatrices();
        void set(const Eigen::MatrixXd& newPositions);
        void update(const Eigen::VectorXd&,  const unsigned int&);
        void reset(const unsigned int&);
        void acceptState(const unsigned int&);

    protected:
        std::vector<double(Quantumdot::*)(const unsigned int&, const unsigned
                int&, const unsigned int&)> variationalDerivativeFunctionList;

        void checkIfFullShell();
};

#endif /* QUANTUMDOT_H */

#endif
