#pragma once
#include "Math/field.h"

namespace LatticeIO {

    class InputConf{
    public:
        static void readConf(GluonField& lattice, int confNum);
        static void readConf(GluonField& lattice, const char *inputFile);
        static void readSubLattice(GluonField& lattice, int confNum);
        static void setInputDir(std::string inputDir);

        static void getInputList(std::vector<std::string>& inputConfList);
    private:
        static int m_linkSize;
        static std::string m_inputDir;
    };

} // end LatticeIO
