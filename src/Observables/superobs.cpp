#include "Observables/superobs.h"
#include "Observables/observable.h"
#include "Math/su3.h"
#include "ParallelTools/parallel.h"
#include "Math/lattice.h"
#include <cstdio>
#include <cmath>

int leviCivit(int i, int j, int k, int l){
    return ( (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12 );
}

// CONSTRUCTOR
SuperObs::SuperObs() : Observable(){
    m_name = "SuperObs";
}

// INITIALIZER
void SuperObs::initObservable(GluonField* field){
    m_field = field;
    trace.allocate(m_field->m_size);
}

void SuperObs::superGatherResults(){
    MPI_Allreduce(&plaq, &plaq, 1, MPI_DOUBLE, MPI_SUM, Parallel::cartCoordComm());
    MPI_Allreduce(&energy, &energy, 1, MPI_DOUBLE, MPI_SUM, Parallel::cartCoordComm());
    MPI_Allreduce(&topc, &topc, 1, MPI_DOUBLE, MPI_SUM, Parallel::cartCoordComm());
}

// COMPUTE THE VALUE OF THE PLAQUETTE ON THE WHOLE LATTICE
void SuperObs::compute(){
    m_value = 0;
    plaqSum.setSU3Zero();
    topc = 0;
    energy = 0;
    // loop through the sublattice

    Lattice<SU3> shiftmu, shiftnu, tempmm, tempmp, temppm;

    for(int mu = 0; mu < 4; mu++){
    for(int nu = mu+1; nu < 4; nu++){

        shiftmu = shift((*m_field)[nu], mu, 1);
        shiftnu = shift((*m_field)[mu], nu, 1);


        clovers[mu][nu]  = (*m_field)[mu];
        clovers[mu][nu] *= shiftmu;
        clovers[mu][nu] *= adj(shiftnu);
        clovers[mu][nu] *=  adj((*m_field)[nu]);


        tempmm   = adj(shiftnu);
        tempmm  *= adj((*m_field)[nu]);
        tempmm  *= (*m_field)[mu];
        tempmm  *= shiftmu;;


        tempmp   = shiftnu;
        tempmp  *= adj(shiftmu);
        tempmp  *= adj((*m_field)[mu]);
        tempmp  *= (*m_field)[nu];


        temppm   = adj((*m_field)[mu]);
        temppm  *= (*m_field)[nu];
        temppm  *= shiftnu;
        temppm  *= adj(shiftmu);


        plaqSum += sum(clovers[mu][nu]);
        clovers[mu][nu] += shift(shift(tempmm, mu, -1), nu, -1);
        clovers[mu][nu] -= shift(tempmp, nu, -1);
        clovers[mu][nu] -= shift(temppm, mu, -1);


/*
        // compute G_mu_nu
        clovers[mu][nu]  =  (*m_field)[mu]
                         *  shift((*m_field)[nu], mu, 1)
                         *  adj(shift((*m_field)[mu], nu, 1))
                         *  adj((*m_field)[nu]);

        plaqSum += sum(clovers[mu][nu]);
        clovers[mu][nu] -=  adj(shift((*m_field)[mu], mu, -1))
                         *  shift((*m_field)[nu], mu, -1)
                         *  shift(shift((*m_field)[mu],nu,1),mu,-1)
                         *  adj((*m_field)[nu]);

        clovers[mu][nu] +=  adj(shift((*m_field)[mu], mu, -1))
                         *  adj(shift(shift((*m_field)[nu],nu,-1),mu,-1))
                         *  shift(shift((*m_field)[mu],nu,-1),mu,-1)
                         *  shift((*m_field)[nu], nu, -1);

        clovers[mu][nu] -=   (*m_field)[mu]
                         *  adj(shift(shift((*m_field)[nu],nu,-1),mu,1))
                         *  adj(shift((*m_field)[mu], nu, -1))
                         *  shift((*m_field)[nu], nu, -1);
*/

        clovers[mu][nu] -= adj(clovers[mu][nu]);
        for(int site = 0; site < clovers[mu][nu].sites; site++){
            double tr = clovers[mu][nu].at(site).imagTrace()/3.0;
            for(int i = 1; i < 18; i+=8)
                clovers[mu][nu].at(site).mat[i] -= tr;
        }


        energy += sum(realTrace((clovers[mu][nu]*clovers[mu][nu])));// + sum(realTrace((adj(clovers[nu][mu])*adj(clovers[nu][mu]))));

   }}


    int mu = 0;
    for(int nu = 1; nu < 4; nu++){
        int rho = nu%3;
        rho++;
        int sig = rho%3;
        sig++;
        if(rho < sig)
            topc -= sum(realTrace((clovers[mu][nu]*clovers[rho][sig])));
        else
            topc -= sum(realTrace((clovers[mu][nu]*adj(clovers[sig][rho]))));
    }
    plaq = plaqSum.realTrace()  / 18.0 / m_field->m_size[0] / m_field->m_size[1] /  m_field->m_size[2] /  m_field->m_size[3] / Parallel::activeProcs();
    topc /= ( 16*16*M_PI*M_PI );
    superGatherResults();
}


