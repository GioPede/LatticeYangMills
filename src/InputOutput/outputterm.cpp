#include <cstdio>
#include "InputOutput/outputterm.h"
#include "Math/lattice.h"
#include "ParallelTools/parallel.h"
#include "Observables/observable.h"

namespace LatticeIO {

/*!
 *   \brief prints basic information about the system to terminal
 */ 
    void OutputTerm::printInitialConditions(){
        if(Parallel::rank() == 0){
            printf("Total Lattice Size: %2i %2i %2i %2i\n", Parallel::latticeFullSize()[0],Parallel::latticeFullSize()[1],
                                                            Parallel::latticeFullSize()[2], Parallel::latticeFullSize()[3]);
            printf("Sublattice Size:    %2i %2i %2i %2i\n\n", Parallel::latticeSubSize()[0],Parallel::latticeSubSize()[1],
                                                              Parallel::latticeSubSize()[2],Parallel::latticeSubSize()[3]);
            printf("Number of Processors: \t\t%i\n", Parallel::numProcs());
            printf("Number of Active Processors: \t%i\n\n", Parallel::activeProcs());
        }
    }

/*!
 *   \brief print information about the system during thermalization
 *   \param step the curren thermalization step
 *   \param obsList a vector of Observables that will be computed and written to file
 *   \param acceptRatio the current acceptance ratio of the Metropolis algorithm
 */ 
    void OutputTerm::printThermStep(int step, std::vector<Observable*>& obsList, double acceptRatio){
        if(Parallel::rank() == 0){
            printf("Thermalization Step %4i\n", step);
            for(auto& obs : obsList)
                printf("\t%s\t = %lf \n", obs->getName(), obs->value());
            printf("\tAccept. Ratio = %lf\n", acceptRatio);
        }
    }

/*!
 *   \brief print information about the system during the gauge field generation
 *   \param confNum the number of the current configuration that has been generated 
 *   \param obsList a vector of Observables that will be computed and written to file
 *   \param acceptRatio the current acceptance ratio of the Metropolis algorithm
 */ 
    void OutputTerm::printGenStep(int confNum, std::vector<Observable*>& obsList, double acceptRatio){
        if(Parallel::rank() == 0){
            printf("Configuration Number %4i\n", confNum);
            for(auto& obs : obsList)
                printf("\t%s\t = %lf \n", obs->getName(), obs->value());
            printf("\tAccept. Ratio = %lf\n", acceptRatio);
        }
    }

/*!
 *   \brief print information about the system of a given configuration
 *   \param confNum the number of the current configuration
 *   \param obsList a vector of Observables that will be computed and written to file
 */ 
    void OutputTerm::writeObs(int confNum, std::vector<Observable*>& obsList){
        if(Parallel::rank() == 0){
            printf("Configuration Number %4i\n", confNum);
            //or(auto& obs : obsList)
            //    printf("\t%s\t = %lf \n", obs->getName(), obs->value());
            printf("\t%s\t = %g \n", "Plaq.", obsList[0]->plaq);
            printf("\t%s\t = %g \n", "Topc.", obsList[0]->topc);
            printf("\t%s\t = %g \n", "Energ.", obsList[0]->energy);
            }

    }

/*!
 *   \brief print information about the system during the gradient flow steps
 *   \param flowTime the current flow time of the integration
 *   \param obsList a vector of Observables that will be computed and written to file
 */ 
    void OutputTerm::writeFlowObs(double flowTime, std::vector<Observable*>& obsList){
        if(Parallel::rank() == 0){
            printf("Flow Time %2.2lf\n", flowTime);
            //for(auto& obs : obsList)
            //    printf("\t%s\t = %lf \n", obs->getName(), obs->value());
            printf("\t%s\t = %g \n", "Plaq.", obsList[0]->plaq);
            printf("\t%s\t = %g \n", "Topc.", obsList[0]->topc);
            printf("\t%s\t = %g \n", "Energ.", obsList[0]->energy);
        }
    }

} // end LatticeIO
