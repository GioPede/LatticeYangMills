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

/*! \file      gaugefieldfactory.h 
 *  \brief     Contains the definition of the GaugeFieldFactory App derived 
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

/*! \class     GaugeFieldFactory
 *  \brief     Implementation of the GaugeFieldFactory App class.
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 * 
 *  The GaugeFieldFactory class implements the Metropolis algorithm to 
 *  generate an ensemble of gauge field configurations. It requires and Action
 *  to be specified and the length of the Markov-Chain to simulate.
 */
class GaugeFieldFactory : public App {
public:
    GaugeFieldFactory(int MCSteps, int thermSteps,
                      int NConf, double epsilon,
                      std::string startType);
    // main functions
    void generateConfigurations();

    // setters-getters
    std::vector<double>& getObsValues() { return m_obsValues; }
    void execute();

private:
    // member classes
    std::vector<double> m_obsValues;

    // Private methods
    void initialize();

    // Input variables
    int m_MCSteps;
    int m_correlationSteps;
    int m_thermSteps;
    int m_confs;
    std::string m_startType;
    double m_epsilon;


    long int m_accepted = 0;
    long int m_updates = 0;
    long m_seed;
    const char* m_outDir;
    void MCUpdate();
    void updateLink(int x, int y, int z, int t, int mu);
    SU3 newLink;
    void computeObservables();
    void thermalize();
    void sampleConf();
    void thermalizeTime();
    void sampleConfTime();
};
