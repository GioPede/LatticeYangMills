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

/*! \file      wilsonflow.cpp
 *  \brief     Contains the implementation of the WilsonFlow class 
 *             methods
 *  \author    Giovanni Pederiva
 *  \version   1.0
 *  \date      2017-2018
 *  \copyright MIT License.
 */

#include "Utils/clusterspecifier.h"
#ifndef LACONIA
#include <mpi/mpi.h>
#else
#include <mpi.h>
#endif
#include <vector>
#include <cstdio>
#include <string>
#include <ctime>
#include <cmath>
#include <random>
#include "lqcd.h"


/*!
 *   \brief Executes the App. First initializes the system then applies the 
 *   gradient flow to the configurations
 */  
void WilsonFlow::execute(){
    initialize();
    flowConfigurations();
}

/*!
 *   \brief constructor of the class
 *   \param tauFinal the final flow time to be reached by the integration
 *   \param epsilon the integration step size
 */  
WilsonFlow::WilsonFlow(double tauFinal, double epsilon){
    m_epsilon = epsilon;
    m_tauFinal = tauFinal;
}

/*!
 *   \brief runs the initializers of the subclasses.
 */  
void WilsonFlow::initialize(){
    for(auto& obs : m_obs)
        obs->initObservable(m_lat);
    m_act->initAction(m_lat);
    m_obsValues.resize(m_obs.size());

    LatticeIO::InputConf::getInputList(m_inputConfList);
    LatticeIO::OutputObs::initialize(m_obs);
    LatticeIO::OutputTerm::printInitialConditions();
}

/*!
 *   \brief creates the GluonField object and the temporary lattice
 *   for the action derivative.
 */  
void WilsonFlow::createLattice(std::array<int,4> latticeSize){
    m_size = latticeSize;
    m_lat = new GluonField(m_size);
    m_Z = new GluonField(m_size);
}

/*!
 *   \brief applies the gradient flow to all the configurations given as input
 */  
void WilsonFlow::flowConfigurations(){
    // check that current processor should be active
    int confNum = 0;
    if(Parallel::isActive()){
        for(auto& conf : m_inputConfList){
            LatticeIO::InputConf::readConf(*m_lat, conf.c_str());
            MPI_Barrier(Parallel::cartCoordComm());
            computeObservables();
            LatticeIO::OutputTerm::writeObs(confNum, m_obs);
            applyWilsonFlow(confNum);
            confNum++;
        }
    }
}

/*!
 *   \brief applies the gradient flow to one configuration
 *   \param confNum the numerical label of the configuration to flow
 */  
void WilsonFlow::applyWilsonFlow(int confNum){
    std::vector<std::vector<double>> flowObsMatrix;
    flowObsMatrix.resize(int(m_tauFinal/m_epsilon));
    for(auto& slice : flowObsMatrix)
        slice.resize(4);
    flowObsMatrix[0][0] = 0.0;
    flowObsMatrix[0][1] =  m_obs[0]->plaq;
    flowObsMatrix[0][2] = m_obs[0]->topc;
    flowObsMatrix[0][3] = m_obs[0]->energy;
    for(int t = 1; t < flowObsMatrix.size(); t++){
        flowStep();
        flowObsMatrix[t][0] = flowObsMatrix[t-1][0] + m_epsilon;
        computeObservables();
        //for(int i = 0; i < m_obs.size(); i++)
        //    flowObsMatrix[t][i+1] = m_obs[i]->value();
        flowObsMatrix[t][1] =  m_obs[0]->plaq;
        flowObsMatrix[t][2] = m_obs[0]->topc;
        flowObsMatrix[t][3] = m_obs[0]->energy;
        LatticeIO::OutputTerm::writeFlowObs(flowObsMatrix[t][0], m_obs);
    }
    LatticeIO::OutputObs::writeFlowObs(confNum, m_obs, flowObsMatrix);
}

/*!
 *   \brief applies one step of the integration scheme, described in paper
 */  
void WilsonFlow::flowStep(){
    // compute Z0
    for(int mu = 0; mu < 4; mu++)
        (*m_Z)[mu] = m_act->computeDerivative(mu);

    // compute W1
    for(int mu = 0; mu < 4; mu++){
        (*m_lat)[mu] = exp((*m_Z)[mu] * m_epsilon * (1.0/4.0))* (*m_lat)[mu];
    }

    // compute Z1
    for(int mu = 0; mu < 4; mu++)
        (*m_Z)[mu] = m_act->computeDerivative(mu)*(8.0/9.0)
                   - (*m_Z)[mu]*(17.0/36.0);


    // compute W2
    for(int mu = 0; mu < 4; mu++){
        (*m_lat)[mu] = exp((*m_Z)[mu]*m_epsilon)* (*m_lat)[mu];
    }

    // compute Z2
    for(int mu = 0; mu < 4; mu++)
        (*m_Z)[mu] = m_act->computeDerivative(mu)*(3.0/4.0)
                   - (*m_Z)[mu];


    // compute V_t+eps
    for(int mu = 0; mu < 4; mu++){
        (*m_lat)[mu] = exp((*m_Z)[mu] *m_epsilon)* (*m_lat)[mu];
    }
}


/*!
 *   \brief computes all observables in the Observable vector
 */  
void WilsonFlow::computeObservables(){
    for(auto& obs : m_obs)
        obs->compute();
}