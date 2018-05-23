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

/*! \file      wilsonflow.h 
 *  \brief     Contains the definition of the WilsonFlow App derived 
 *             class
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#pragma once
#include <vector>
#include <random>
#include "Apps/app.h"

/*! \class     WilsonFlow
 *  \brief     Implementation of the WilsonFlow App class.
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  The WilsonFlow class implements Gradiend Flow method as suggested by 
 *  M.Lüscher in the article https://arxiv.org/abs/1006.4518 
 *  The integration scheme is a third order Runge-Kutta method for Lie
 *  algebras, so called Runge-Kutta Munthe-Kaas RK-MK method.
 */

class WilsonFlow : public App {
    public:
        WilsonFlow(double tauFinal, double epsilon);
        // main functions
        void flowConfigurations();

        // setters-getters
        void setAction(Action* action);
        void addObservable(Observable* observable);
        //class Point& getLatticeSite(int x, int y, int z, int t);
        std::array<int, 4>& getSize() { return m_size; }
        std::vector<double>& getObsValues() { return m_obsValues; }
        std::vector<Observable*>& getObs() { return m_obs; }
        void createLattice(std::array<int,4> latticeSize);
        void execute();
        void initialize();


    private:
        // flow lattices
        GluonField* m_Z = nullptr;
        std::vector<double> m_obsValues;
        std::vector<std::string> m_inputConfList;

        SU3 omega;
        // member variables
        double m_epsilon;
        double m_tauFinal;
        void computeObservables();
        void flowStep();
        void applyWilsonFlow(int confNum);
};
