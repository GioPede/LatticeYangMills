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

/*! \file      inputconf.h 
 *  \brief     Contains classes for reading lattices from binary files
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once

#include "Math/field.h"

namespace LatticeIO {

    /*! \class     InputConf
     *  \brief     Class for reding lattices from binary files
     *  \author    Giovanni Pederiva
     *  \version   1.0
     *  \date      2017-2018
     *  \copyright MIT License.
     * 
     *  This class mainly serves the purpose of providing an interface 
     *  for reading GluonField objects from binary format.
     */
    class InputConf{
    public:
        /// Reads a binary file from a int marker for the name
        static void readConf(GluonField& lattice, int confNum);
        /// Reads a binary file given its path into a lattice object
        static void readConf(GluonField& lattice, const char *inputFile);
        /// Reads one binary file for each processor, assuming it is split
        static void readSubLattice(GluonField& lattice, int confNum);
        /// Sets the input directory
        static void setInputDir(std::string inputDir);
        /// Returns a list of all valid input files in a the input directory
        static void getInputList(std::vector<std::string>& inputConfList);
    
    private:
        /// Contains the size in bytes of a 4 links on a lattice site 
        static int m_linkSize;
        /// The path of the input directory
        static std::string m_inputDir;
    };

} // end LatticeIO
