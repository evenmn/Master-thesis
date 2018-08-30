#ifndef TESTWAVEBASIS_H
#define TESTWAVEBASIS_H

#include <Eigen/Dense>

class testWaveBasis {
    private:
        unsigned int m_dim;

    public:
        testWaveBasis ();
        virtual ~testWaveBasis ();

       void setup();
};

#endif /* TESTWAVEBASIS_H */
