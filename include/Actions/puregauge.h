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

/*! \file      puregauge.h 
 *  \brief     Contains the definition of the PureGauge action derived class
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once

#include "Actions/action.h"
#include "Math/latticemath.h"

//! Initializer for the PureGauge Action class.
/*!
 *   beta The $/beta$ value of the action
 */
class PureGauge : public Action {
public:
    PureGauge(GluonField* lattice, double beta) {
              m_lat = lattice; m_beta = beta; }
    PureGauge(double beta);
    double compute(int x, int y, int z, int t, int mu, SU3& newLink);
    void   computeStaples(int mu);
    Lattice<SU3>    computeDerivative(int mu);
    void   computeStaplez(GluonField* lattice);

    void computeOtherStaples(int x, int y, int z, int t, int mu);
private:
    double m_beta;
    GluonField* m_staplesLattice = nullptr;

    // auxiliary variables
    GluonField* m_constPart;
    Lattice<SU3> m_omega;
    Lattice<SU3> trace;
};

