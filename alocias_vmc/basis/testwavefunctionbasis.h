#ifndef TESTWAVEFUNCTIONBASIS_H
#define TESTWAVEFUNCTIONBASIS_H

#include <Eigen/Dense>

class TestWavefunctionBasis {
    private:
        unsigned int m_dim;

    public:
        TestWavefunctionBasis ();
        virtual ~TestWavefunctionBasis ();

       void setup();
};

#endif /* TESTWAVEFUNCTIONBASIS_H */
