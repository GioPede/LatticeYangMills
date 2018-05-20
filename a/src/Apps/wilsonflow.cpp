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


void WilsonFlow::execute(){
    initialize();
    flowConfigurations();
}

WilsonFlow::WilsonFlow(double tauFinal, double epsilon){
    m_epsilon = epsilon;
    m_tauFinal = tauFinal;
}

void WilsonFlow::initialize(){
    for(auto& obs : m_obs)
        obs->initObservable(m_lat);
    m_act->initAction(m_lat);
    m_obsValues.resize(m_obs.size());

    LatticeIO::InputConf::getInputList(m_inputConfList);
    LatticeIO::OutputObs::initialize(m_obs);
    LatticeIO::OutputTerm::printInitialConditions();
}


void WilsonFlow::createLattice(std::array<int,4> latticeSize){
    m_size = latticeSize;
    m_lat = new GluonField(m_size);
    m_Z = new GluonField(m_size);
}

// MAIN FUNCTION OF CLASS. GENERATES GAUGE FIELD CONFIGURATION USING METROPILIS'
// ALGORITHM AND SAVES THEM FOR FUTURE ANALISYS.
void WilsonFlow::flowConfigurations(){
    // check that current processor should be active
    int confNum = 0;
    double epsilon = m_epsilon;
    if(Parallel::isActive()){
        for(auto& conf : m_inputConfList){
            LatticeIO::InputConf::readConf(*m_lat, conf.c_str());
	    MPI_Barrier(Parallel::cartCoordComm());
		if(Parallel::rank() == 0) printf("ready: size of array %i\n\n", (*m_lat)[3].lattice.size());
	    computeObservables();
            LatticeIO::OutputTerm::writeObs(confNum, m_obs);
            applyWilsonFlow(confNum, epsilon);
            confNum++;
        }
    }
}

void WilsonFlow::applyWilsonFlow(int confNum, double epsilon){
    std::vector<std::vector<double>> flowObsMatrix;
    flowObsMatrix.resize(int(m_tauFinal/epsilon));
    for(auto& slice : flowObsMatrix)
        slice.resize(4);
    flowObsMatrix[0][0] = 0.0;
    //for(int i = 0; i < m_obs.size(); i++)
    //    flowObsMatrix[0][i+1] = m_obs[i]->value();
    flowObsMatrix[0][1] =  m_obs[0]->plaq;
    flowObsMatrix[0][2] = m_obs[0]->topc;
    flowObsMatrix[0][3] = m_obs[0]->energy;
    for(int t = 1; t < flowObsMatrix.size(); t++){
        flowStep(epsilon);
        flowObsMatrix[t][0] = flowObsMatrix[t-1][0] + epsilon;
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

void WilsonFlow::flowStep(double epsilon){
    // compute Z0
    for(int mu = 0; mu < 4; mu++)
        (*m_Z)[mu] = m_act->computeDerivative(mu);

    // compute W1
    for(int mu = 0; mu < 4; mu++){
        (*m_lat)[mu] = exp((*m_Z)[mu] * epsilon * (1.0/4.0))* (*m_lat)[mu];
    }

    // compute Z1
    for(int mu = 0; mu < 4; mu++)
        (*m_Z)[mu] = m_act->computeDerivative(mu)*(8.0/9.0)
                   - (*m_Z)[mu]*(17.0/36.0);


    // compute W2
    for(int mu = 0; mu < 4; mu++){
        (*m_lat)[mu] = exp((*m_Z)[mu]*epsilon)* (*m_lat)[mu];
    }

    // compute Z2
    for(int mu = 0; mu < 4; mu++)
        (*m_Z)[mu] = m_act->computeDerivative(mu)*(3.0/4.0)
                   - (*m_Z)[mu];


    // compute V_t+eps
    for(int mu = 0; mu < 4; mu++){
        (*m_lat)[mu] = exp((*m_Z)[mu] *epsilon)* (*m_lat)[mu];
    }
}



void WilsonFlow::computeObservables(){
    for(auto& obs : m_obs)
        obs->compute();
}

/*
// GETTERS AND SETTERS
Point& WilsonFlow::getLatticeSite(int x, int y, int z, int t){
    return (*m_lat)(x,y,z,t);
}*/

void WilsonFlow::setAction(Action* action){
    m_act = action;
}

void WilsonFlow::addObservable(Observable *observable){
    m_obs.push_back(observable);
    m_obsValues.push_back(0.0);
}
