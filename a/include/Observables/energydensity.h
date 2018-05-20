#pragma once
#include <vector>
#include <string>
#include "observable.h"
#include "Math/latticemath.h"

class EnergyDensity : public Observable {
public:
    EnergyDensity();
    void initObservable(GluonField *lattice);
    void compute();

private:
    std::array<int, 4> m_size;
    // auxiliary variables
    double m_norm;
    SU3 Gmn;
};
