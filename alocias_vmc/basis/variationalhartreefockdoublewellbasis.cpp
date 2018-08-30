#include "variationalhartreefockdoublewellbasis.h"

VariationalHartreeFockDoubleWellBasis::VariationalHartreeFockDoubleWellBasis()
    : Cartesian(), DWC() {
    /* default constructor */
} // end constructor

VariationalHartreeFockDoubleWellBasis::~VariationalHartreeFockDoubleWellBasis()
{
} // end deconstructor

void VariationalHartreeFockDoubleWellBasis::setup(unsigned int dim) {
    /* initiate states */
    DWC::setup(dim);

    Cartesian::setup(2*DWC::C.rows(), dim);
    Cartesian::restructureStates();
} // end function setup 
