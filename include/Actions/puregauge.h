#pragma once
#include "Actions/action.h"
#include "Math/latticemath.h"

// Implementation of the Pure Gauge Wilson Action
class PureGauge : public Action {
public:
    PureGauge(GluonField* lattice, double beta) {
              m_lat = lattice; m_beta = beta; }
    PureGauge(double beta);
    double compute(int x, int y, int z, int t, int mu, SU3& newLink);
    void   computeStaples(int mu);
    Lattice<SU3>    computeDerivative(int mu);
    void   computeStaplez(GluonField* lattice);

    void computeOtherStaples(int x, int y, int z, int t, int mu);
private:
    double m_beta;
    GluonField* m_staplesLattice = nullptr;

    // auxiliary variables
    GluonField* m_constPart;
    Lattice<SU3> m_omega;
    Lattice<SU3> trace;
};

