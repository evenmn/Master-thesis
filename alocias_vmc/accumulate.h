#ifndef ACCUMULATE_H
#define ACCUMULATE_H

#include <eigen3/Eigen/Dense>

struct Accumulative {
    /* struct containing accumulative values used in VMC. Overloads for /=, *=,
     * =, +=, -=, +, - and * are given. */
    double acceptance;
    double energy;
    double energySquared;
    double potentialEnergy;
    double kineticEnergy;
    Eigen::VectorXd psiDerivative;
    Eigen::VectorXd psiDerivativeTimesEnergy;
    Eigen::VectorXd psiDerivativeTimesEnergySquared;

    Accumulative& operator << (const size_t& n) {
        /* initialize vector with size n */
        this->acceptance = 0;
        this->energy = 0;
        this->energySquared = 0;
        this->potentialEnergy = 0;
        this->kineticEnergy = 0;
        this->psiDerivative = Eigen::VectorXd::Zero(n);
        this->psiDerivativeTimesEnergy = Eigen::VectorXd::Zero(n);
        this->psiDerivativeTimesEnergySquared = Eigen::VectorXd::Zero(n);

        return *this;
    } // end operator << overload

    template<typename T> Accumulative& operator /= (const T& value) {
        /* overload divide operator to work as if struct Accumulative
         * was a vector */
        this->acceptance /= value;
        this->energy /= value;
        this->energySquared /= value;
        this->potentialEnergy /= value;
        this->kineticEnergy /= value;
        this->psiDerivative /= value;
        this->psiDerivativeTimesEnergy /= value;
        this->psiDerivativeTimesEnergySquared /= value;

        return *this;
    } // end operator /= overload
    
    template<typename T> Accumulative& operator *= (const T& value) {
        /* overload multiplication operator to work as if struct
         * Accumulative was a vector */
        this->acceptance *= value;
        this->energy *= value;
        this->energySquared *= value;
        this->potentialEnergy *= value;
        this->kineticEnergy *= value;
        this->psiDerivative *= value;
        this->psiDerivativeTimesEnergy *= value;
        this->psiDerivativeTimesEnergySquared *= value;

        return *this;
    } // end operator *= overload
    
    template<typename T> Accumulative& operator = (const T& value) {
        /* overload = to set all elements in struct to value */
        this->acceptance = value;
        this->energy = value;
        this->energySquared = value;
        this->potentialEnergy = value;
        this->kineticEnergy = value;
        Eigen::VectorXd tmpVec =
            Eigen::VectorXd::Constant(this->psiDerivative.size(), value);
        this->psiDerivative = tmpVec;
        this->psiDerivativeTimesEnergy = tmpVec;
        this->psiDerivativeTimesEnergySquared = tmpVec;

        return *this;
    } // end operator = overload
    
    template<typename T> Accumulative& operator += (const T& value) {
        /* overload += to add value to all elements in struct */
        this->acceptance += value;
        this->energy += value;
        this->energySquared += value;
        this->potentialEnergy += value;
        this->kineticEnergy += value;
        this->psiDerivative = this->psiDerivative.array() + value;
        this->psiDerivativeTimesEnergy = this->psiDerivativeTimesEnergy.array()
            + value;
        this->psiDerivativeTimesEnergySquared =
            this->psiDerivativeTimesEnergySquared.array() + value;

        return *this;
    } // end operator += overload
    
    template<typename T> Accumulative& operator -= (const T& value) {
        /* overload -= to subtract value from all elements in struct */
        this->acceptance -= value;
        this->energy -= value;
        this->energySquared -= value;
        this->potentialEnergy -= value;
        this->kineticEnergy -= value;
        this->psiDerivative = this->psiDerivative.array() - value;
        this->psiDerivativeTimesEnergy = this->psiDerivativeTimesEnergy.array()
            - value;
        this->psiDerivativeTimesEnergySquared =
            this->psiDerivativeTimesEnergySquared.array() - value;

        return *this;
    } // end operator -= overload

    Accumulative& operator + (Accumulative& obj) {
        /* overload + to add corresponding values in obj to
         * m_accumulative */
        this->acceptance += obj.acceptance;
        this->energy += obj.energy;
        this->energySquared += obj.energySquared;
        this->potentialEnergy += obj.potentialEnergy;
        this->kineticEnergy += obj.kineticEnergy;
        this->psiDerivative += obj.psiDerivative;
        this->psiDerivativeTimesEnergy += obj.psiDerivativeTimesEnergy;
        this->psiDerivativeTimesEnergySquared +=
            obj.psiDerivativeTimesEnergySquared;

        return *this;
    } // end operator + overload
    
    Accumulative& operator - (Accumulative& obj) {
        /* overload - to subtract corresponding values in obj to
         * m_accumulative */
        this->acceptance -= obj.acceptance;
        this->energy -= obj.energy;
        this->energySquared -= obj.energySquared;
        this->potentialEnergy -= obj.potentialEnergy;
        this->kineticEnergy -= obj.kineticEnergy;
        this->psiDerivative -= obj.psiDerivative;
        this->psiDerivativeTimesEnergy -= obj.psiDerivativeTimesEnergy;
        this->psiDerivativeTimesEnergySquared -=
            obj.psiDerivativeTimesEnergySquared;

        return *this;
    } // end operator - overload
    
    Accumulative& operator * (Accumulative& obj) {
        /* overload multiplication operator to work as if struct
         * Accumulative was a vector */
        this->acceptance *= obj.acceptance;
        this->energy *= obj.energy;
        this->energySquared *= obj.energySquared;
        this->potentialEnergy *= obj.potentialEnergy;
        this->kineticEnergy *= obj.kineticEnergy;
        this->psiDerivative *= obj.psiDerivative;
        this->psiDerivativeTimesEnergy *= obj.psiDerivativeTimesEnergy;
        this->psiDerivativeTimesEnergySquared *=
            obj.psiDerivativeTimesEnergySquared;

        return *this;
    } // end operator * overload
};

#endif /* ACCUMULATE_H */
