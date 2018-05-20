#include <cstdio>
#include <cmath>
#include "Actions/puregauge.h"
#include "Math/latticemath.h"
#include "Utils/latticeunits.h"
#include "ParallelTools/parallel.h"
#include "Utils/clusterspecifier.h"
#ifndef LACONIA
#include <mpi/mpi.h>
#else
#include <mpi.h>
#endif

// COMPUTE THE ACTION DIFFERENCE ON SINGLE LINK AFTER UPDATE, TO BE USED IN
// METROPOLIS ALGORITHM
PureGauge::PureGauge(double beta){
    m_beta = beta;
    m_staplesLattice = new GluonField(Parallel::latticeSubSize());
    m_constPart = new GluonField(Parallel::latticeSubSize());
    m_omega.allocate(Parallel::latticeSubSize());
    trace.allocate(Parallel::latticeSubSize());
    LatticeUnits::initialize(m_beta);
}

double PureGauge::compute(int x, int y, int z, int t, int mu, SU3& newLink){
    //return -m_beta / 3.0 * getMultSumTrace(newLink - (*m_lat)(x,y,z,t)[mu], m_constPart);

}


// COMPUTE THE SUM OF THE STAPLES OF A GIVEN LINK
void PureGauge::computeStaples(int mu){
    setToZero(m_omega);
    Lattice<SU3> shiftmu, tempStaple;
    for(int nu = 0; nu < 4; nu++){
        if(nu != mu){
            shiftmu = shift((*m_lat)[nu],mu,1);
            tempStaple  = shiftmu;
            tempStaple *= adj( shift((*m_lat)[mu],nu,1));
            tempStaple *= adj( (*m_lat)[nu]);
            m_omega += tempStaple;

            tempStaple  = adj(shiftmu);
            tempStaple *= adj( (*m_lat)[mu]);
            tempStaple *= (*m_lat)[nu];

            m_omega += shift(tempStaple, nu, -1);
        }
    }
}


Lattice<SU3> PureGauge::computeDerivative(int mu){
    computeStaples(mu);
    m_omega = adj(m_omega) * (adj((*m_lat)[mu]));
    m_omega = (adj(m_omega)-m_omega);

    for(int site = 0; site < m_omega.sites; site++){
        double tr = m_omega.at(site).imagTrace()/3.0;
        for(int i = 1; i < 18; i+=8)
            m_omega.at(site).mat[i] -= tr;

        for(int i = 0; i < 18; i+=2){
            double temp = m_omega.at(site).mat[i];
            m_omega.at(site).mat[i] = -m_omega.at(site).mat[i+1]*0.5;
            m_omega.at(site).mat[i+1] = temp*0.5;
        }
    }
    return m_omega;
}
