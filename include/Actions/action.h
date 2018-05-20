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

/*! \file      action.h 
 *  \brief     Contains the definition of the Action prototype
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once

#include "Math/field.h"

template <typename T> class Lattice;
struct SU3;

/*! \class     Action
 *  \brief     Prototype for the Action class group
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  The Action class is one of the key elements of the program. It specifies
 *  the acceptance of a Metropolis test and the evolution according to flow 
 *  equations. It is of virtual type, so an implementation is always required.
 *  It mainly offers the possibility of retrieving the action around a single 
 *  link or the derivative of the action around all links.
 */
class Action{
public:
    /// Computes the action difference around a given link for a suggested move
    virtual double compute(int x, int y, int z, int t, int mu, SU3& newLink) = 0;
    /// Computes the sum of the staples around all links
    virtual void   computeStaples(int mu) = 0;
    /// Computes the derivative in a given direction around all links
    virtual Lattice<SU3>    computeDerivative(int mu) = 0;
    /// Links the action to a given field and initializes
    void initAction(GluonField* field);

protected:
    /// The GluonField object the action is linked to
    GluonField* m_lat = nullptr;
};
