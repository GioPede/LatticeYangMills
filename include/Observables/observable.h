#pragma once
#include <string>
#include "Math/field.h"

// Prototype for class Observable. Every implementation of a physical observable
// should inherit from this class.
class Observable{
public:
    // virtual methods
    virtual void compute() = 0;
    virtual void initObservable(GluonField *field) = 0;
    const char* getName();
    double value();

    double plaq, energy, topc;
protected:
    void gatherResults();
    GluonField* m_field = nullptr;
    double m_value;
    std::string m_name;
};

