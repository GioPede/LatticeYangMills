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

/*! \file      superobs.h 
 *  \brief     Contains the definition of the SuperObs class
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once
#include <vector>
#include <string>
#include "Observables/observable.h"
#include "Math/latticemath.h"

/*! \class     SuperObs
 *  \brief     Class to compute the plaquette, the energy density and the 
 *             topological charge
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  Currently used observable class. Computes the clover for the Field Strength 
 *  Tensor and uses it to compute the plaquette, the energy density and the 
 *  topological charge. This should be changed...
 */
class SuperObs : public Observable {
public:
    SuperObs();
    void initObservable(GluonField *field);
    void compute();

private:
    Lattice<SU3> trace;
    // auxiliary variables
    double m_norm;

    std::array<std::array<Lattice<SU3>,4>,4> clovers;
    SU3 plaqSum;

    std::array<int,24> mu;
    std::array<int,24> nu;
    std::array<int,24> rho;
    std::array<int,24> sig;
    void superGatherResults();
};
