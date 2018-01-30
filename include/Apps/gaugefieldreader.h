#pragma once
#include <vector>
#include <array>
#include "Apps/app.h"

class GaugeFieldReader : public App {
public:
    GaugeFieldReader() {}
    // main functions
    void initGFR();
    void sampleConfigurations();

    // setters-getters
    void addObservable(Observable* observable);
    const char* getOutDir() { return m_outDir; }
    //class Point& getLatticeSite(int x, int y, int z, int t);
    std::array<int, 4>& getSize() { return m_size; }
    std::vector<double>& getObsValues() { return m_obsValues; }
    std::vector<Observable*>& getObs() { return m_obs; }
    void execute();


private:
    // member classes
    GluonField* m_lat = nullptr;
    std::vector<Observable*> m_obs;
    std::vector<double> m_obsValues;
    std::vector<std::string> m_inputConfList;
    std::array<int, 4> m_size;


    // member variables
    const char* m_outDir;
};
