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

/*! \file      outputobs.h 
 *  \brief     Contains classes for output to file of observables values
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 *
 *  This file defines the class needed to file the values of the observables at
 *  different steps in the execution. It handles through MPI the acces to files
 *  depending on rank.
 */

#pragma once

#include <vector>
#include <string>

#include "Observables/observable.h"

namespace LatticeIO {

    /*! \class     OutputObs
     *  \brief     Class for output to file of observables values
     *  \author    Giovanni Pederiva
     *  \version   1.0
     *  \date      2017-2018
     *  \copyright MIT License.
     * 
     *  This file defines the class needed to file the values of the 
     *  observables at different steps in the execution. It handles through 
     *  MPI the acces to files depending on rank. Provides output for 
     *  generation and flow steps.
     */
    class OutputObs{
    public:
        /// Opens the files and writes needed headers
        static void initialize(std::vector<Observable*>& obsList);
        /// Writes the list of observables to file at given generation point
        static void writeObs(std::vector<Observable*>& obsList, int MCSteps);
        /// Writes the observables for given flow step
        static void writeFlowObs(int confNum, std::vector<Observable*>& obsList, std::vector<std::vector<double>>& obsMatrix);
        /// Sets the output directory for the file output
        static void setOutputDir(std::string outputDir);

    private:
        /// The file object for target output
        static FILE* m_file;
        /// The name of the output file
        static char fileName[1024];
        /// The output directory
        static std::string m_outputDir;

    };

} // end LatticeIO
