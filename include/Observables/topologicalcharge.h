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

/*! \file      topologicalcharge.h 
 *  \brief     Contains the definition of the TopologicalCharge observable
 *  \author    Giovanni Pederiva
 *  \version   0.1
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once
#include <vector>
#include <string>
#include "Observables/observable.h"
#include "Math/latticemath.h"

/*! \class     TopologicalCharge
 *  \brief     Implementation of the Topological Charge operator class
 *  \author    Giovanni Pederiva
 *  \version   0.1
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  Old implementation of the Topological Charge operator class using the 
 *  clover definition of the Field Stength Tensor
 */
class TopologicalCharge : public Observable {
public:
    TopologicalCharge();
    void initObservable(Lattice *lattice);
    void compute();

private:
    std::array<int, 4> m_size;
    // auxiliary variables
    double m_norm;
    SU3 Gmn, Grs;
};
