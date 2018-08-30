#ifndef SLATERJASTROW_H
#define SLATERJASTROW_H

#include <eigen3/Eigen/Dense>

#ifdef EXPNQS
    #include "jastrow/expnqs.h"
    class ExpNQS;
    using Jastrow = ExpNQS;
#endif

#ifdef PADENQS
    #include "jastrow/padenqs.h"
    class PadeNQS;
    using Jastrow = PadeNQS;
#endif

#ifdef PADEJASTROW 
    #include "jastrow/padejastrow.h"
    class PadeJastrow;
    using Jastrow = PadeJastrow;
#endif

#ifdef RBMJASTROW
    #include "jastrow/rbmjastrow.h"
    class RBMJastrow;
    using Jastrow = RBMJastrow;
#endif

#include "slater.h"

class SlaterJastrow : public Slater, public Jastrow {
    private:
        Eigen::MatrixXd m_jastrowGradientMatrix, m_oldJastrowGradientMatrix;

        unsigned int m_numVariationalParameters;

        // expand macro which checks if reset and accept functions exist in
        // jastrow class (at compile time) and compiles in call specified in
        // .cpp
        MEM_FUNC(resetGradient, has_resetGradient, jastrowResetGradient);
        MEM_FUNC(acceptGradient, has_acceptGradient, jastrowAcceptGradient);
        MEM_FUNC(calculateGradient, has_calculateGradient,
                jastrowCalculateGradient);

    public:
        SlaterJastrow (const unsigned int&, const unsigned int&, const
                Eigen::VectorXd&);
        virtual ~SlaterJastrow ();

        double wavefunctionRatio();
        double wavefunctionRatio(const unsigned int&);
        void calculateGradient();
        void calculateGradient(const unsigned int&);

        double kineticEnergy();

        void setVariationalDerivatives();

        void initializeMatrices();
       
        void acceptState(const unsigned int&);
        void reset(const unsigned int&);
        void update(const Eigen::VectorXd&, const unsigned int&);
        void resetGradient(const unsigned int&);
        void acceptGradient(const unsigned int&);
        void set(const Eigen::MatrixXd);
        
        const double &getJastrowGradient(const unsigned int&, const unsigned
                int&) const;
        const Eigen::Ref<const Eigen::RowVectorXd> getJastrowGradient(const
                unsigned int&) const;
        Eigen::RowVectorXd getGradient(const unsigned int&);
        const Eigen::MatrixXd &getJastrowGradientMatrix() const;
        const double &getOldJastrowGradient(const unsigned int&, const unsigned
                int&) const;
        const Eigen::Ref<const Eigen::RowVectorXd> getOldJastrowGradient(const
                unsigned int&) const;
        const Eigen::MatrixXd &getOldJastrowGradientMatrix() const;
};

#endif /* SLATERJASTROW_H */
