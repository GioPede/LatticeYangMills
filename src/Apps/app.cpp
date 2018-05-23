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

/*! \file      app.cpp
 *  \brief     Contains the implementation of the App prototype class methods
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#include "lqcd.h"


//! Initializes the GluonField object given a sublattice size
/*!
  \param latticeSize 4-component int array for the size of the sublattice
*/
void App::createLattice(std::array<int,4> latticeSize){
    m_size = latticeSize;
    m_lat = new GluonField(m_size);
}

//! Initializes the Action object pointer given an instance of it
/*!
  \param action a pointer to an Action  instance
*/
void App::setAction(Action* action){
    m_act = action;
}

//! Adds a Observable class pointer to the vector of observables of the App
/*!
  \param observable a pointer to an Observable instance
*/
void App::addObservable(Observable *observable){
    m_obs.push_back(observable);
}
