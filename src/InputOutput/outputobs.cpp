#include <cstdio>
#include <boost/filesystem.hpp>
#include "Observables/observable.h"
#include "InputOutput/outputobs.h"
#include "Math/lattice.h"
#include "ParallelTools/parallel.h"

namespace LatticeIO {

    char OutputObs::fileName[1024];
    FILE* OutputObs::m_file;
    std::string OutputObs::m_outputDir;


/*!
 *   \brief sets the output directory static variable
 *   \param outputDir a string for the path of the output directory
 */ 
    void OutputObs::setOutputDir(std::string outputDir){
        m_outputDir = outputDir;
        boost::filesystem::create_directory(m_outputDir);
        boost::filesystem::create_directory(m_outputDir+"/flow");
    }


/*!
 *   \brief initializes the output file for the observables
 *   \param obsList a vector of Observables that will be computed and written to file
 */ 
    void OutputObs::initialize(std::vector<Observable*>& obsList){
        if(Parallel::rank() == 0){
            sprintf(fileName, "%s/observables.dat", m_outputDir.c_str());
            m_file = fopen(fileName, "w+");
            fprintf(m_file, "Step\t");
            for(auto& obs : obsList)
                fprintf(m_file, "%s\t", obs->getName());
            fprintf(m_file, "\n");
            fclose(m_file);
        }
    }

/*!
 *   \brief writes the current state of the observables to file, given a numerical tag
 *   \param obsList a vector of Observables that will be computed and written to file
 *   \param MCSteps the tag of the current Monte Carlo step
 */ 
    void OutputObs::writeObs(std::vector<Observable*>& obsList, int MCSteps){
        if(Parallel::rank() == 0){
            m_file = fopen(fileName, "a");
            fprintf(m_file, "%i\t", MCSteps);
            for(auto& obs : obsList)
                fprintf(m_file, "%lf\t", obs->value());
            fprintf(m_file, "\n");
            fclose(m_file);
        }
    }

/*!
 *   \brief writes the matrix of observables for each flowtime to file, given a numerical tag
 *   \param confNum the tag of the current configuration that's being flowed
 *   \param obsList a vector of Observables that will be computed and written to file
 *   \param obsMatrix a matrix containing all the observables computed at every flow time
 */ 
    void OutputObs::writeFlowObs(int confNum, std::vector<Observable*>& obsList, std::vector<std::vector<double>>& obsMatrix){
        if(Parallel::rank() == 0){
            char flowFileName [1024];
            sprintf(flowFileName, "%s/flow/conf%04d.dat", m_outputDir.c_str(), confNum);
            FILE* flowOut = fopen(flowFileName, "w+");
            fprintf(flowOut, "Flow Time\t");
            //for(auto& obs : obsList)
            //    fprintf(flowOut, "%s\t", obs->getName());
            fprintf(flowOut, "%s\t", "Plaq.");
            fprintf(flowOut, "%s\t", "Topc.");
            fprintf(flowOut, "%s\t", "Energ.");
            fprintf(flowOut, "\n");

            for(auto& flowStep : obsMatrix){
                for(auto& obs : flowStep)
                    fprintf(flowOut, "%lf\t", obs);
                fprintf(flowOut, "\n");
            }
            fclose(flowOut);
        }
    }

} // end LatticeIO
