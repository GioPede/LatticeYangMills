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

/*! \file      outputterm.h 
 *  \brief     Contains classes for output to standard command line interface
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 *
 *  This file defines the class needed for ouputting to standard command line 
 *  interface, handling the multiprocess nature of the code. Namely it only 
 *  allows one processor to output any data, preventing the others from doing 
 *  anything.
 */

#pragma once

#include <vector>

#include "Math/lattice.h"
#include "Observables/observable.h"

namespace LatticeIO {

    /*! \class     OutputTerm
     *  \brief     Class for standard output
     *  \author    Giovanni Pederiva
     *  \version   1.0
     *  \date      2017-2018
     *  \copyright MIT License.
     * 
     *  This class handles the standard terminal output dealing with the
     *  multi processor nature of the program. It checks on the MPI rank in 
     *  order to decide whether to output or not. It also provides basic 
     *  tools to print information about the system as it evolves.
     */
    class OutputTerm{
    public:
        /// Prints basic info about the system
        static void printInitialConditions();
        /// Prints the status of the current object list for a given conf
        static void writeObs(int confNum, std::vector<Observable*>& obsList);
        /// Prints information about a thermalization step
        static void printThermStep(int step, std::vector<Observable*>& obsList, double acceptRatio);
        /// Prints information about a generation step
        static void printGenStep(int confNum, std::vector<Observable*>& obsList, double acceptRatio);
        /// Prints information about a flow time evolution step
        static void writeFlowObs(double flowTime, std::vector<Observable*>& obsList);
    };

} // end LatticeIO
