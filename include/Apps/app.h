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

/*! \file      app.h 
 *  \brief     Contains the definition of the App prototype
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once
#include <array>
#include <vector>
#include "Math/field.h"
#include "Observables/observable.h"

template <typename T> class Lattice;
class Action;
struct SU3;


/*! \class     App
 *  \brief     Prototype for the App class group
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  The App class represents the main user interface of the program. It is
 *  initialized by the JsonInput class when an input file is passed. The app
 *  holds the Action class, the main GluonField object and the Observable 
 *  vector. 
 *  The generic execute() method is used to start the execution of the main 
 *  algorithm.
 */
class App{
public:
    void setAction(Action *action);
    virtual void createLattice(std::array<int,4> latticeSize);
    void addObservable(Observable* observable);

    virtual void execute() = 0;

protected:
    /// Pointer to the Action instance 
    Action* m_act = nullptr;
    /// Pointer to the GluonField instance 
    GluonField* m_lat = nullptr;
    /// Size of the sublattice
    std::array<int,4> m_size;
    /// Vector of pointers to Observable instances
    std::vector<Observable*> m_obs;

};
