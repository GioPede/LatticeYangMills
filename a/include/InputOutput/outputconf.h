#pragma once
#include "Math/field.h"

namespace LatticeIO {

    class OutputConf{
    public:
        static void writeConf(GluonField& lattice, int confNum);
        static void writeSubLattice(GluonField& lattice, int confNum);
        static void setOutputDir(std::string outputDir);
    private:
        static int m_linkSize;
        static std::string m_outputDir;

    };

} // end LatticeIO
