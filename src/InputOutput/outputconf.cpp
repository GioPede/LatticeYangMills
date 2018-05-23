/******************************************************************************
* 
* MIT License
* 
* Copyright (c) 2018 Giovanni Pederiva
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in 
* all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
******************************************************************************/

/*! \file      outputconf.cpp
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#include <boost/filesystem.hpp>
#include "Utils/clusterspecifier.h"
#ifndef LACONIA
#include <mpi/mpi.h>
#else
#include <mpi.h>
#endif
#include <cstdio>

#include "InputOutput/outputconf.h"
#include "Math/lattice.h"
#include "ParallelTools/parallel.h"


namespace LatticeIO {

    
    int OutputConf::m_linkSize = 72 * sizeof(double);
    std::string OutputConf::m_outputDir;

/*!
 *   \brief sets the output directory static variable
 *   \param outputDir a string for the path of the output directory
 */ 
    void OutputConf::setOutputDir(std::string outputDir){
        m_outputDir = outputDir;
        boost::filesystem::create_directory(m_outputDir);
    }

    // // Write a single file for every sublattice (for testing)
    // void OutputConf::writeSubLattice(GluonField& lattice, int confNum){
    //     char fileName [1024];
    //     sprintf(fileName, "%s/r%03du%03d", m_outputDir.c_str(), Parallel::rank(), confNum);
    //     FILE* output = fopen(fileName, "wb");
    //     for(int t = 0; t < Parallel::latticeSubSize()[3]; t++){
    //     for(int z = 0; z < Parallel::latticeSubSize()[2]; z++){
    //     for(int y = 0; y < Parallel::latticeSubSize()[1]; y++){
    //     for(int x = 0; x < Parallel::latticeSubSize()[0]; x++){
    //         //fwrite(lattice(x,y,z,t).m_links, sizeof(double), 72, output);
    //     }}}}
    //     fclose(output);
    // }

/*!
 *   \brief writes a lattice configuration in a given field object and a numerical tag
 *   \param lattice the GluonField to write in the default output folder
 *   \param confNum the numerical tag for the configuration to read
 */     void OutputConf::writeConf(GluonField& lattice, int confNum){
        MPI_File output;
        MPI_Offset startPointT, startPointZ, startPointY, startPointX;
        int volumeX = 1;
        int volumeY = Parallel::latticeFullSize()[0] * volumeX;
        int volumeZ = Parallel::latticeFullSize()[1] * volumeY;
        int volumeT = Parallel::latticeFullSize()[2] * volumeZ;

        // Create filename and open
        char fileName [1024];
        sprintf(fileName, "%s/conf%04d.bin", m_outputDir.c_str(), confNum);
        Parallel::openFile(output, fileName);

        for(int t = 0; t < Parallel::latticeSubSize()[3]; t++){
            startPointT = volumeT * ( Parallel::rankCoord()[3]*Parallel::latticeSubSize()[3] + t);
            for(int z = 0; z < Parallel::latticeSubSize()[2]; z++){
            startPointZ = startPointT + volumeZ * ( Parallel::rankCoord()[2]*Parallel::latticeSubSize()[2] + z);
                for(int y = 0; y < Parallel::latticeSubSize()[1]; y++){
                startPointY = startPointZ + volumeY * ( Parallel::rankCoord()[1]*Parallel::latticeSubSize()[1] + y);
                    for(int x = 0; x < Parallel::latticeSubSize()[0]; x++){
                        startPointX = startPointY + volumeX * ( Parallel::rankCoord()[0]*Parallel::latticeSubSize()[0] + x);
                        startPointX *= m_linkSize;
                        //MPI_File_write_at(output, startPointX, lattice(x,y,z,t).m_links, 72, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
        Parallel::closeFile(output);
    }

} // end LatticeIO
