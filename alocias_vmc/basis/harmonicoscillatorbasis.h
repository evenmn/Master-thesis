#ifndef QUANTUMDOTBASIS_H
#define QUANTUMDOTBASIS_H

#include <Eigen/Dense>

#include "../methods.h"
#include "../basisfunctions/cartesian.h"

class QuantumdotBasis : public Cartesian {
    private:
        unsigned int m_dim, E;
        void printDescription();
    
    public:
        QuantumdotBasis (const unsigned int, const unsigned int);
        QuantumdotBasis ();
        virtual ~QuantumdotBasis ();
        
        void setup(const unsigned int, const unsigned int);
        void printStates();
};

#endif /* QUANTUMDOTBASIS_H */
