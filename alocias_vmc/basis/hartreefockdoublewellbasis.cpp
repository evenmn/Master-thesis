#include "hartreefockdoublewellbasis.h"

HartreeFockDoubleWellBasis::HartreeFockDoubleWellBasis() : Cartesian(), DWC() {
    /* default constructor */
} // end constructor

HartreeFockDoubleWellBasis::~HartreeFockDoubleWellBasis() {
} // end deconstructor

void HartreeFockDoubleWellBasis::setup(unsigned int dim) {
    /* initiate states */
    DWC::setup(dim);
    
    Cartesian::setup(2*DWC::C.rows(), dim);
    Cartesian::restructureStates();
} // end function setup 
