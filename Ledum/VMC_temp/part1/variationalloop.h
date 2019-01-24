#pragma once
#include "variationalmc.h"


class VariationalLoop {
    private:
        VariationalMC vmc;

        int         cycles              = 0;
        double      energy              = 0;
        double      alpha               = 0;
        double      beta                = 0;
        double      molecularDistance   = 0;
        double      accepted            = 0;
        arma::vec   energyVarGrad       = zeros<vec>(2);


    public:
        VariationalLoop();
        void initialize_processes(int);
        void run();
        void setAlphaBeta(double alpha, double beta);
        void setMolecularDistance(double R);
        void setNumberOfCycles(int n);
        void setNumberOfMonteCarloCycles(int n);


        double getEnergy() const;
        void setEnergy(double value);
        double getAlpha() const;
        void setAlpha(double value);
        double getBeta() const;
        void setBeta(double value);
        double getAccepted() const;
        void setAccepted(double value);
        double getEnergyVarGrad(int i) const;
        void setEnergyVarGrad(const arma::vec& value);
};

