#ifndef HARTREEFOCKDOUBLEWELLBASIS_H
#define HARTREEFOCKDOUBLEWELLBASIS_H

#include "../basisfunctions/cartesian.h"
#include "../basisfunctions/dwc.h"

#include <Eigen/Dense>

class HartreeFockDoubleWellBasis : public Cartesian, public DWC {
    private:
        unsigned int m_dim;

    public:
        HartreeFockDoubleWellBasis ();
        virtual ~HartreeFockDoubleWellBasis ();

        void setup(unsigned int);
};

#endif /* HARTREEFOCKDOUBLEWELLBASIS_H */
