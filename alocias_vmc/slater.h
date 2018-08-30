#ifndef SLATER_H
#define SLATER_H

#ifdef HARMONICOSCILLATOR 
    #include "wavefunctions/harmonicoscillator.h"
    class Quantumdot;
    using WF = Quantumdot;
#endif

#ifdef HARTREEFOCKDOUBLEWELL
    #include "wavefunctions/hartreefockdoublewell.h"
    class HartreeFockDoubleWell;
    using WF = HartreeFockDoubleWell;
#endif

#ifdef VARIATIONALHARTREEFOCKDOUBLEWELL
    #include "wavefunctions/variationalhartreefockdoublewell.h"
    class VariationalHartreeFockDoubleWell;
    using WF = VariationalHartreeFockDoubleWell;
#endif

#ifdef HARTREEFOCK 
    #include "wavefunctions/hartreefock.h"
    class HartreeFock;
    using WF = HartreeFock;
#endif

#ifdef VARIATIONALHARTREEFOCK
    #include "wavefunctions/variationalhartreefock.h"
    class VariationalHartreeFock;
    using WF = VariationalHartreeFock;
#endif

#ifdef TESTWAVE 
    #include "wavefunctions/testwave.h"
    class testWave;
    using WF = testWave;
#endif

#include <eigen3/Eigen/Dense>
#include <iostream>

#include "hasmemfunc.h"


class Slater : public WF {
    private:
        // expand above macro for functions set, reSetAll, initializeMatrices,
        // update, reset, resetGradient, acceptState and acceptGradient. The
        // last name in each call to MEM_FUNC is the corresponding function
        // used to call the respective function (if they exist) in the source
        // file.
        MEM_FUNC(set, has_set, wfset);
        MEM_FUNC(reSetAll, has_reSetAll, wfreSetAll);
        MEM_FUNC(initializeMatrices, has_initializeMatrices,
                wfinitializeMatrices);
        MEM_FUNC(update, has_update, wfupdate);
        MEM_FUNC(reset, has_reset, wfreset);
        MEM_FUNC(resetGradient, has_resetGradient, wfresetGradient);
        MEM_FUNC(acceptState, has_acceptState, wfacceptState);
        MEM_FUNC(acceptGradient, has_acceptGradient, wfacceptGradient);

        void setIndices(const unsigned int&);

        void setDistances();
        void setDistances(const unsigned int&);
        void setDistances(const unsigned int&, const unsigned int&);

        void resetDistance(const unsigned int&);

        void setGradient(const unsigned int&, const unsigned int&);

    protected:
        unsigned int invIdx, waveIdx, spanIdx;

        Eigen::MatrixXd m_newPositions, m_oldPositions, m_newDistances,
            m_oldDistances;

        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
            m_newDistanceMatrix, m_oldDistanceMatrix;

        Eigen::MatrixXd m_gradient, m_oldGradient;

        Eigen::MatrixXd m_newWavefunctionMatrix, m_oldWavefunctionMatrix,
            m_newWavefunctionMatrixInverse, m_oldWavefunctionMatrixInverse;
        
        double determinantRatio;


    public:
        Slater (const unsigned int&, const unsigned int&, const
                Eigen::VectorXd&);
        virtual ~Slater ();
        
        unsigned int m_dim, m_numParticles, halfSize;

        Eigen::VectorXd m_parameters, m_firstDerivativesParameters;

        void initializeMatrices();

        void setSize();
        void setSize(const unsigned int&);

        double wavefunctionRatio();
        double wavefunctionRatio(const double&);

        void calculateWavefunction(const unsigned int&);
        
        void setWavefunction();
        void setWavefunction(const unsigned int&);

        void setVariationalDerivatives();
        
        void calculateGradient();
        void calculateGradient(const unsigned int&);

        double laplacian();

        void update(const Eigen::VectorXd&, const unsigned int&);
        void reset(const unsigned int&);

        void updateWavefunction(const unsigned int&);
        void updateWavefunction();
        
        void resetGradient(const unsigned int&);
        
        void set(const Eigen::MatrixXd&);
        void reSetAll();
        
        void acceptState(const unsigned int&);
        
        void acceptGradient(const unsigned int&);

        double kineticEnergy();
       
        const unsigned int& getSpan() const;
        const unsigned int &getDimension() const;
        const unsigned int &getNumberOfParticles() const;
        const Eigen::VectorXd& getParameters() const;
        const Eigen::VectorXd& getVariationalDerivatives() const;
        const double &getNewPosition(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getNewPosition(const unsigned
                int&) const;
        const Eigen::MatrixXd &getNewPositions() const;
        const double &getOldPosition(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getOldPosition(const unsigned
                int&) const;
        const Eigen::MatrixXd &getOldPositions() const;
        const double &getGradient(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getGradient(const unsigned
                int&) const;
        const Eigen::Ref<const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
              1>> &getGradient3D(const unsigned int&) const;
        const Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
            &getGradient3DMatrix() const;
        const Eigen::MatrixXd &getGradientMatrix() const;
        const Eigen::Ref<const Eigen::MatrixXd> getInverse() const;
        const Eigen::Ref<const Eigen::MatrixXd> getOldInverse() const;
        const Eigen::MatrixXd &getInverseMatrix() const;
        const Eigen::MatrixXd &getOldInverseMatrix() const;
        const Eigen::MatrixXd &getWavefunctionMatrix() const;
        const Eigen::MatrixXd &getOldWavefunctionMatrix() const;
        const Eigen::Ref<const Eigen::RowVectorXd> getWavefunction(const
                unsigned int& i) const;
        const Eigen::Ref<const Eigen::RowVectorXd> getOldWavefunction(const
                unsigned int& i) const;
        const double &getWavefunction(const unsigned int&, const unsigned int&)
            const;
        const double &getOldWavefunction(const unsigned int&, const unsigned
                int&) const;
        const double &getNewDistance(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getNewDistance(const unsigned
                int&) const;
        const Eigen::MatrixXd &getNewDistanceMatrix() const;
        const double &getOldDistance(const unsigned int&, const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::RowVectorXd> getOldDistance(const unsigned
                int&) const;
        const Eigen::MatrixXd &getOldDistanceMatrix() const;
        const Eigen::Ref<const Eigen::RowVectorXd> getNewDistanceVector(const
                unsigned int&, const unsigned int&);
        const Eigen::Ref<const Eigen::RowVectorXd> getOldDistanceVector(const
                unsigned int&, const unsigned int&);
};

typedef Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>
EigenMatVecXd;

#endif /* SLATER_H */
