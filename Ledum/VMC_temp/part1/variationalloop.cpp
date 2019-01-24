#include <iostream>
#include <fstream>
#include <iomanip>
#include "variationalloop.h"


double VariationalLoop::getEnergy() const
{
    return energy;
}

void VariationalLoop::setEnergy(double value)
{
    energy = value;
}

double VariationalLoop::getAlpha() const
{
    return alpha;
}

void VariationalLoop::setAlpha(double value)
{
    alpha = value;
}

double VariationalLoop::getBeta() const
{
    return beta;
}

void VariationalLoop::setBeta(double value)
{
    beta = value;
}

double VariationalLoop::getAccepted() const
{
    return accepted;
}

void VariationalLoop::setAccepted(double value)
{
    accepted = value;
}

double VariationalLoop::getEnergyVarGrad(int i) const
{
    return energyVarGrad(i);
}

void VariationalLoop::setEnergyVarGrad(const arma::vec& value)
{
    energyVarGrad = value;
}

VariationalLoop::VariationalLoop() {

}


void VariationalLoop::run() {
    double startClock, finishClock, energy, accepted;
    vec    energyVarGrad = zeros<vec>(4);
    vec    varGrad       = zeros<vec>(2);
    vec    alphaBeta     = zeros<vec>(2);


    alphaBeta(0) = this->alpha;
    alphaBeta(1) = this->beta;

    for (int i = 1; i < this->cycles+1; i++) {
        startClock = clock();

        // Run one metropolis loop.
        energyVarGrad = vmc.runMetropolis(alphaBeta(0), alphaBeta(1));

        energy        = energyVarGrad(0);
        varGrad(0)    = energyVarGrad(1);
        varGrad(1)    = energyVarGrad(2);
        accepted      = energyVarGrad(3);

        if (accepted > 0.95) {
            // Compute new values of alpha / beta.
            alphaBeta     += (1.0/i) * varGrad * energy;
        } else {
            //cout << "\n\n\n\n\n\naccepted  = "<< accepted << " < 0.95 \n\n\n\n\n";
        }


        if (alphaBeta(0) < 0) {
            alphaBeta(0) = 0.01;
        } if (alphaBeta(1) < 0) {
            alphaBeta(1) = 0.01;
        }
        this->alpha = alphaBeta(0);
        this->beta = alphaBeta(1);
        this->energy = energy;
        this->energyVarGrad(0) = varGrad(0);
        this->energyVarGrad(1) = varGrad(1);
        this->accepted = accepted;

        finishClock = clock();
        //cout << "   * Time usage       = " << (finishClock - startClock) / 1000000.0 << " [s] "<< endl;
        //cout << "   * Alpha            = " << alphaBeta(0) << endl;
        //cout << "   * Beta             = " << alphaBeta(1) << endl;
        //cout << "VarGrad=" << varGrad << endl;
        //cout << endl << endl;
    }
}

void VariationalLoop::setAlphaBeta(double alpha, double beta) {
    this->alpha = alpha;
    this->beta = beta;
}

void VariationalLoop::setMolecularDistance(double R) {
    this->molecularDistance = R;
    this->vmc.setMolecularDistance(R);
}

void VariationalLoop::setNumberOfCycles(int n) {
    this->cycles = n;
}

void VariationalLoop::setNumberOfMonteCarloCycles(int n) {
    this->vmc.setNumberOfMonteCarloCycles(n);
}

