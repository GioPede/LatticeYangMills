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

/*! \file      parallel.h 
 *  \brief     Utilities for parallelization
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once
#include "Utils/clusterspecifier.h"
#ifndef LACONIA
#include <mpi/mpi.h>
#else
#include <mpi.h>
#endif
#include <array>

/*! \class     Parallel
 *  \brief     Class to manage the parallelization scheme
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  Utility class to manage the parallel architecture of the program
 *  and interface each process with the overall job. Handles local
 *  coordinates, nearest neighbors and communicators
 */


class Parallel{
// Public interface of the Parallel class
public:
    // Init-Final methods
    static void initialize();
    static void createGeometry(std::array<int, 4> latticeSize,
                               std::array<int, 4> subLatticeSize);
    static void finalize();

    // Data members
    static int rank() { return m_rank; }
    static int numProcs()  { return m_numProcs; }
    static int activeProcs() { return m_activeProcs; }
    static bool isActive() { return m_isActive; }
    static MPI_Comm cartCoordComm() { return m_cartCoordComm; }
    static std::array<int, 4>& subBlocks() { return m_subBlocks; }
    static std::array<int, 4>& rankCoord() { return m_rankCoord; }
    static std::array<int, 4>& latticeSubSize() { return m_latticeSubSize; }
    static std::array<int, 4>& latticeFullSize() { return m_latticeFullSize; }
    static std::array<int, 4>& parity(){ return m_parity; }

    // Neighbor getters
    static int getNeighbor(int direction, int sign);
    static int getSecondNeighbor(int direction1, int sign1,
                                 int direction2, int sign2);

    // File IO
    static void openFile (MPI_File& file, const char *fileName);
    static void closeFile(MPI_File& file);

private:
    // Members made private for safety
    static int m_rank;
    static int m_numProcs;
    static int m_activeProcs;
    static bool m_isActive;
    static MPI_Comm m_cartCoordComm;
    static std::array<int, 4> m_subBlocks;
    static std::array<int, 4> m_rankCoord;
    static std::array<int, 4> m_latticeSubSize;
    static std::array<int, 4> m_latticeFullSize;
    static std::array<int, 4> m_parity;
    static MPI_Comm m_cartCoords;
    static std::array< std::array<int, 2>, 4> m_neighbor;
    static std::array< std::array<std::array<std::array<int, 2>, 4>, 2>, 4> m_secondNeighbor;

    // Init-Only methods
    static void createNeighborLists();
    static void assignNeighbor(int direction);
    static void assignSecondNeighbor(int dir1, int dir2);
    static void MyMPI_Cart_shift2(MPI_Comm comm, int dir1, int shift1,
                                                 int dir2, int shift2,
                                                 int &source, int &dest);
};
