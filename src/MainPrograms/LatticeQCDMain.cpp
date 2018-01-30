#include <cstdio>
#include <chrono>
#include "lqcd.h"

int main(int argn, char* argv[]) {

    // Start timer
    auto start = std::chrono::system_clock::now();

    // Initialize MPI
    Parallel::initialize();

    // Read and validate the input file
    LatticeIO::JsonInput::parseInput(argn, argv);

    // Create Application
    App* lqcd;
    LatticeIO::JsonInput::initApp(lqcd);

    // Execute app
    lqcd->execute();

    // Finalize MPI and delete objects
    Parallel::finalize();

    // Stop timer and output
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime = end-start;
    if(Parallel::rank() == 0)
        printf("Total time taken: %.2fs\n", elapsedTime.count());
/*
    Parallel::initialize();
    Parallel::createGeometry(std::array<int,4> {8,8,8,8}, std::array<int,4> {4,4,4,4});

    Lattice field(std::array<int,4> {4,4,4,4});

    LatticeIO::InputConf::setInputDir("/home/giovanni/Desktop/Lattice/conf/8_16");
    LatticeIO::InputConf::readConf(field, "/home/giovanni/Desktop/Lattice/conf/8_16/conf0002.bin");

    Observable* obs = new SuperObs();
    obs->initObservable(&field);

    auto start = std::chrono::system_clock::now();
    for(int i = 0; i < 1000; i++){
        obs->compute();
    }
    auto end = std::chrono::system_clock::now();
    if(Parallel::rank() == 0){
        printf("Plaq = %g\n", obs->plaq);
        printf("Ener = %g\n", obs->energy);
        printf("Topc = %g\n", obs->topc);
    }

    // Finalize MPI and delete objects
    Parallel::finalize();

    // Stop timer and output
    //auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime = end-start;
    if(Parallel::rank() == 0)
        printf("Total time taken: %.2fs\n", elapsedTime.count());
*/
    return 0;
}
