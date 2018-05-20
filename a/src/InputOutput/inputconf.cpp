#include "Utils/clusterspecifier.h"
#ifndef LACONIA
#include <mpi/mpi.h>
#else
#include <mpi.h>
#endif
#include <cstdio>
#include <algorithm>
#include<boost/filesystem.hpp>
#include "InputOutput/inputconf.h"
#include "Math/lattice.h"
#include "ParallelTools/parallel.h"
#include "Utils/clusterspecifier.h"

namespace LatticeIO {

    double ReverseFloat( const double inFloat )
    {
       double retVal;
       char *floatToConvert = ( char* ) & inFloat;
       char *returnFloat = ( char* ) & retVal;

       // swap the bytes into a temporary buffer
       returnFloat[0] = floatToConvert[7];
       returnFloat[1] = floatToConvert[6];
       returnFloat[2] = floatToConvert[5];
       returnFloat[3] = floatToConvert[4];
       returnFloat[4] = floatToConvert[3];
       returnFloat[5] = floatToConvert[2];
       returnFloat[6] = floatToConvert[1];
       returnFloat[7] = floatToConvert[0];
       return retVal;
    }

    double ReverseDouble( const double inDouble ){
       double retVal;
       char *doubleToConvert = ( char* ) & inDouble;
       char *returnDouble = ( char* ) & retVal;

       // swap the bytes into a temporary buffer
       returnDouble[0] = doubleToConvert[7];
       returnDouble[1] = doubleToConvert[6];
       returnDouble[2] = doubleToConvert[5];
       returnDouble[3] = doubleToConvert[4];
       returnDouble[4] = doubleToConvert[3];
       returnDouble[5] = doubleToConvert[2];
       returnDouble[6] = doubleToConvert[1];
       returnDouble[7] = doubleToConvert[0];
       return retVal;

    }


    int InputConf::m_linkSize = 18 * sizeof(double);
    std::string InputConf::m_inputDir;

    void InputConf::setInputDir(std::string inputDir){
        m_inputDir = inputDir;
    }

    //if (boost::filesystem::is_directory(inputPath))

    void InputConf::getInputList(std::vector<std::string>& inputConfList){
        boost::filesystem::path inputPath = m_inputDir;
        for (boost::filesystem::directory_entry& x : boost::filesystem::directory_iterator(inputPath))
            if(x.path().extension() == ".bin")
                inputConfList.push_back(x.path().string());
		std::sort(inputConfList.begin(), inputConfList.end());
    }

    // Read a single file for every sublattice (for testing)
    void InputConf::readSubLattice(GluonField& lattice, int confNum){
        char fileName [1024];
        sprintf(fileName, "%s/r%03du%03d", OUT_PREFIX, Parallel::rank(), confNum);
        FILE* input = fopen(fileName, "rb");
        for(int t = 0; t < Parallel::latticeSubSize()[3]; t++){
        for(int z = 0; z < Parallel::latticeSubSize()[2]; z++){
        for(int y = 0; y < Parallel::latticeSubSize()[1]; y++){
        for(int x = 0; x < Parallel::latticeSubSize()[0]; x++){
            //fread(lattice(x,y,z,t).m_links, sizeof(double), 72, input);
        }}}}
        fclose(input);
    }

    // Use MPI to read a file containing the lattice
    void InputConf::readConf(GluonField& lattice, int confNum){
        char fileName[1024];
        sprintf(fileName, "%s/conf%04d.bin", OUT_PREFIX, confNum);
        readConf(lattice, fileName);
    }

    // Use MPI to read a file containing the lattice
    void InputConf::readConf(GluonField& lattice, const char* inputFile){
        MPI_File input;
        MPI_Offset startPointT, startPointZ, startPointY, startPointX;
        int volumeX = 1;
        int volumeY = Parallel::latticeFullSize()[0] * volumeX;
        int volumeZ = Parallel::latticeFullSize()[1] * volumeY;
        int volumeT = Parallel::latticeFullSize()[2] * volumeZ;

        // Open file
        Parallel::openFile(input, inputFile);

        for(int t = 0; t < Parallel::latticeSubSize()[3]; t++){
            startPointT = volumeT * ( Parallel::rankCoord()[3]*Parallel::latticeSubSize()[3] + t);
            for(int z = 0; z < Parallel::latticeSubSize()[2]; z++){
            startPointZ = startPointT + volumeZ * ( Parallel::rankCoord()[2]*Parallel::latticeSubSize()[2] + z);
                for(int y = 0; y < Parallel::latticeSubSize()[1]; y++){
                startPointY = startPointZ + volumeY * ( Parallel::rankCoord()[1]*Parallel::latticeSubSize()[1] + y);
                    for(int x = 0; x < Parallel::latticeSubSize()[0]; x++){
                        startPointX = startPointY + volumeX * ( Parallel::rankCoord()[0]*Parallel::latticeSubSize()[0] + x);
                        startPointX *= 4*m_linkSize;

                        // FOR CHROMA
//                        for(int mu = 0; mu < 4; mu++){
//                            for(int i = 0; i < 18; i++){
//                                double a;
//                                MPI_File_read_at(input, startPointX, &a, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
//                                lattice[mu].at(x,y,z,t).mat[i] = ReverseDouble(a);
//                                startPointX += sizeof(double);
//                            }
//                        }

                        // FOR REASONABLE ENDIAN
                        for(int mu = 0; mu < 4; mu++){
                            MPI_File_read_at(input, startPointX, lattice[mu].at(x,y,z,t).mat.data(), 18, MPI_DOUBLE, MPI_STATUS_IGNORE);
                            startPointX += m_linkSize;
                        }
                    }
                }
            }
        }
        Parallel::closeFile(input);
    }

} // end LatticeIO
