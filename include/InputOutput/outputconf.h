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

/*! \file      outputconf.h 
 *  \brief     Contains classes for saving lattices to binary files
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 *
 *  This file mainly serves the purpose of providing an interface 
 *  for saving GluonField objects into binary format.
 */

#pragma once

#include "Math/field.h"

namespace LatticeIO {

    /*! \class     OutputConf
     *  \brief     Class for saving lattices to binary files
     *  \author    Giovanni Pederiva
     *  \version   1.0
     *  \date      2017-2018
     *  \copyright MIT License.
     * 
     *  This object allows saving of GluonField objects as binary files.
     *  It uses MPI routines for handling a single output file with multiple 
     *  processors.
     */
    class OutputConf{
    public:
        /// Saves the given GluonField object to the ouput directory 
        static void writeConf(GluonField& lattice, int confNum);
        /// Saves the GluonField object in a set of files, one for each core
        static void writeSubLattice(GluonField& lattice, int confNum);
        /// Sets the output directory path
        static void setOutputDir(std::string outputDir);
    private:
        /// Contains the size in bytes of a 4 links on a lattice site 
        static int m_linkSize;
        /// The path of the output directory
        static std::string m_outputDir;

    };

} // end LatticeIO
