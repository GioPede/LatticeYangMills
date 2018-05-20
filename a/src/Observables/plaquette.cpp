#include "Observables/plaquette.h"
#include "Observables/observable.h"
#include "Math/su3.h"
#include "ParallelTools/parallel.h"
#include "Math/lattice.h"
#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>

// CONSTRUCTOR
Plaquette::Plaquette() : Observable(){
    m_name = "Plaquette";
}


// INITIALIZER
void Plaquette::initObservable(Lattice* lattice){
    m_lat = lattice;
    m_size = m_lat->getSize();
    m_norm = 1.0 / 18.0 / m_size[0] / m_size[1] / m_size[2] / m_size[3] / Parallel::numProcs();
}

// COMPUTE THE VALUE OF THE PLAQUETTE ON THE WHOLE LATTICE
void Plaquette::compute(){
    // loop through the sublattice
    m_value = 0;
    for(int mu = 0; mu < 4; mu++){
        for(int nu = mu + 1; nu < 4; nu++){
            for(int x = 0; x < m_size[0]; x++){
            for(int y = 0; y < m_size[1]; y++){
            for(int z = 0; z < m_size[2]; z++){
            for(int t = 0; t < m_size[3]; t++){
                m_value += (   (*m_lat)(x,y,z,t)[mu]
                             *   m_lat->shift(x,y,z,t,nu, mu, 1)
                             * ~(m_lat->shift(x,y,z,t,mu, nu, 1))
                             * ~(*m_lat)(x,y,z,t)[nu]
                           ).realTrace();
            }}}}
        }
    }
    m_value *= m_norm;
    gatherResults();
}
