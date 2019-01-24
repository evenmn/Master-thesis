#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <armadillo>
#include "lib.h"
#include "variationalmc.h"
#include "variationalloop.h"

using namespace std;
using namespace arma;


int main(int argc, char* argv[]) {


    int n   = 150;
    vec R   = zeros<vec>(n);
    double Rmin = 3.5;
    double Rmax = 4.5;
    double h    = (Rmax - Rmin) / (n - 1);


    for (int i=0; i < n; i++) {
        R(i) = Rmin + i*h;
    }

    /*
    double  alpha       = 1.0;
    double  beta        = 1.0;
    int     cyclesFirst = 200;
    int     cycles      = 100;
    int     MCCyclesFirst   = (int) 1e8;
    int     MCCycles        = (int) 1e7;
    */
    double  alpha       = 1.0;
    double  beta        = 1.0;
    int     cyclesFirst = 50;
    int     cycles      = 20;
    int     MCCyclesFirst   = (int) 1e7;
    int     MCCycles        = (int) 1e7;

    VariationalLoop loop;
    loop.setAlphaBeta(alpha, beta);

    ofstream outFile;
    outFile.open("VMCData35-45.dat", ios::out);

    for (int i=0; i < n; i++) {

        // Use more cycles for the first run, since our guess for alpha and
        // beta are probably bad.
        if (i==0) {
            loop.setNumberOfCycles(cyclesFirst);
            loop.setNumberOfMonteCarloCycles(MCCyclesFirst);
        } else {
            loop.setNumberOfCycles(cycles);
            loop.setNumberOfMonteCarloCycles(MCCycles);
        }
        loop.setMolecularDistance(R(i));
        loop.run();
        /*                     R    E     a     b    da     db   acc
        sprintf(outString, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f",
                            R(i),
                            loop.getEnergy(),
                            loop.getAlpha(),
                            loop.getBeta(),
                            loop.getEnergyVarGrad(0),
                            loop.getEnergyVarGrad(1),
                            loop.getAccepted());
        outFile << setprecision(15) << outString << endl;*/
        outFile << setprecision(15)
                << R(i)                         << " "
                << loop.getEnergy()             << " "
                << loop.getAlpha()              << " "
                << loop.getBeta()               << " "
                << loop.getEnergyVarGrad(0)     << " "
                << loop.getEnergyVarGrad(1)     << " "
                << loop.getAccepted()           << endl;
        cout    << setprecision(5)
                << R(i)                         << " "
                << loop.getEnergy()             << " "
                << loop.getAlpha()              << " "
                << loop.getBeta()               << " "
                << loop.getEnergyVarGrad(0)     << " "
                << loop.getEnergyVarGrad(1)     << " "
                << loop.getAccepted()           << endl;

    }
    outFile.close();
    return 0;
}


