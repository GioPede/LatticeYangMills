#pragma once
#include <array>
#include <vector>
#include "Math/field.h"
#include "Observables/observable.h"

template <typename T> class Lattice;
class Action;
struct SU3;

class App{
public:
    void setAction(Action *action);
    virtual void createLattice(std::array<int,4> latticeSize);
    void addObservable(Observable* observable);

    // Virtual execution
    virtual void execute() = 0;

protected:
    Action* m_act = nullptr;
    GluonField* m_lat = nullptr;
    std::array<int,4> m_size;
    std::vector<Observable*> m_obs;

};
