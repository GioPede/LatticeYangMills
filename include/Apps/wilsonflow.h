#pragma once
#include <vector>
#include <random>
#include <array>
#include "Apps/app.h"
#include "Math/su3.h"

class WilsonFlow : public App {
    public:
        WilsonFlow(double tauFinal, double epsilon);
        // main functions
        void flowConfigurations();

        // setters-getters
        void setAction(Action* action);
        void addObservable(Observable* observable);
        //class Point& getLatticeSite(int x, int y, int z, int t);
        std::array<int, 4>& getSize() { return m_size; }
        std::vector<double>& getObsValues() { return m_obsValues; }
        std::vector<Observable*>& getObs() { return m_obs; }
        void createLattice(std::array<int,4> latticeSize);
        void execute();
        void initialize();


    private:
        // flow lattices
        GluonField* m_Z = nullptr;
        std::vector<double> m_obsValues;
        std::vector<std::string> m_inputConfList;

        SU3 omega;
        // member variables
        double m_epsilon;
        double m_tauFinal;
        void computeObservables();
        void flowStep(double epsilon);
        void applyWilsonFlow(int confNum, double epsilon);
};
