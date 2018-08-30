#ifndef VARIATIONALHARTREEFOCKDOUBLEWELLBASIS_H
#define VARIATIONALHARTREEFOCKDOUBLEWELLBASIS_H

#include "../basisfunctions/cartesian.h"
#include "../basisfunctions/dwc.h"

#include <Eigen/Dense>

class VariationalHartreeFockDoubleWellBasis : public Cartesian, public DWC {
    private:
        unsigned int m_dim;

    public:
        VariationalHartreeFockDoubleWellBasis ();
        VariationalHartreeFockDoubleWellBasis (unsigned int);
        virtual ~VariationalHartreeFockDoubleWellBasis ();

       void setup(unsigned int);
};

#endif /* VARIATIONALHARTREEFOCKDOUBLEWELLBASIS_H */
