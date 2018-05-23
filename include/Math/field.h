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

/*! \file      field.h 
 *  \brief     Contains the definition of the Field class
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once
#include <vector>
#include <random>
#include "latticemath.h"

/*! \class     Field
 *  \brief     Class that represents the QCD field as an array of 4 SU3 
 *             lattices
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  The Field class is a convenience object that stores 4 lattices of  one 
 *  dattype, one for each euclidean dimension. It allows indexed access to 
 *  the individual components via the operator[].
 */

template <typename T, int N>
class Field{
public:
    // Constructors
    Field(std::array<int, 4> size);
    inline Lattice<T>& operator[](int mu) { return m_lats[mu]; }

    int m_dimensions;
    std::array<int, 4> m_size;
private:
    std::array<Lattice<T>, N> m_lats;
};

template <typename T, int N>
Field<T,N>::Field(std::array<int, 4> size){
    // Allocate memory
    m_size = size;
    m_dimensions = N;
    for(int i = 0; i < N; i++)
        m_lats[i].allocate(m_size);
}

/// Convenience definition of the SU3 case 4 dimensional case of a Field class,
/// the GluonField
typedef Field<SU3,4> GluonField;
