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

/*! \file      gaugefieldreader.h 
 *  \brief     Contains the definition of the GaugeFieldReader App derived 
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

/*! \class     GaugeFieldReader
 *  \brief     Implementation of the GaugeFieldReader App class.
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  The GaugeFieldReader app is a simple program to print individual links of
 *  a given gauge field configuration (an input file).
 */

class GaugeFieldReader : public App {
public:
    GaugeFieldReader() {}
    // main functions
    void initGFR();
    void navigateConf();

    // setters-getters
    void addObservable(Observable* observable);
    const char* getOutDir() { return m_outDir; }
    std::array<int, 4>& getSize() { return m_size; }
    std::vector<double>& getObsValues() { return m_obsValues; }
    std::vector<Observable*>& getObs() { return m_obs; }
    void execute();


private:
    // member classes
    GluonField* m_lat = nullptr;
    std::vector<Observable*> m_obs;
    std::vector<double> m_obsValues;
    std::vector<std::string> m_inputConfList;
    std::array<int, 4> m_size;


    // member variables
    const char* m_outDir;
};
