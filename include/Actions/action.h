#pragma once
#include "Math/field.h"

template <typename T> class Lattice;
struct SU3;

// Prototype for class Action. Every implementation of a physical action should
// inherit from this class.
class Action{
public:
    // virtual methods for the action class
    virtual double compute(int x, int y, int z, int t, int mu, SU3& newLink) = 0;
    virtual void   computeStaples(int mu) = 0;
    virtual Lattice<SU3>    computeDerivative(int mu) = 0;
    virtual void   computeStaplez(GluonField* lattice) = 0;
    void initAction(GluonField* field);

protected:
    GluonField* m_lat = nullptr;
};
