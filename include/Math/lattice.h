#pragma once
#include <vector>
#include <array>
#include <utility>
#include "su3.h"
#include "Utils/clusterspecifier.h"
#ifndef LACONIA
#include <mpi/mpi.h>
#else
#include <mpi.h>
#endif
#include "ParallelTools/parallel.h"

template <typename T>
class Lattice{
public:
    // Default Constructors
    Lattice();
    // Copy Constructors
    Lattice(const Lattice& other) noexcept{
        allocate(other.size);
        for(int i = 0; i < sites; i++)
            lattice[i] = other.lattice[i];
    }
    Lattice(Lattice&& other) noexcept :
        lattice(std::move(other.lattice)),
        sites(std::move(other.sites)),
        size(std::move(other.size)) { }

    // Array-based constructor
    Lattice(std::array<int, 4> sizeArray);
    void allocate(std::array<int, 4> sizeArray);

    // Data
    std::vector<T> lattice;
    std::array<int,4> size;
    int sites;

    // Element accessors
    inline T& at(int x, int y, int z, int t) {
        return lattice[x + size[0]*(y + size[1]*(z + size[2]*t))];
    }

    inline T& at(const std::vector<int>& site){
        return lattice[site[0] + size[0]*(site[1] + size[1]*(site[2] + size[2]*site[3]))];
    }
    inline T& at(const std::array<int, 4>& site){
        return lattice[site[0] + size[0]*(site[1] + size[1]*(site[2] + size[2]*site[3]))];
    }
    inline T& at(int i){
        return lattice[i];
    }
    inline const T& at(int x, int y, int z, int t) const {
        return lattice[x + size[0]*(y + size[1]*(z + size[2]*t))];
    }

    inline const T& at(const std::vector<int>& site) const {
        return lattice[site[0] + size[0]*(site[1] + size[1]*(site[2] + size[2]*site[3]))];
    }
    inline const T& at(const std::array<int, 4>& site) const {
        return lattice[site[0] + size[0]*(site[1] + size[1]*(site[2] + size[2]*site[3]))];
    }
    inline const T& at(int i) const {
        return lattice[i];
    }

    // Copy-assignment
    inline Lattice& operator=(const Lattice& other) noexcept{
        lattice = other.lattice;
        sites = other.sites;
        size = other.size;
        return *this;
    }
    inline Lattice& operator=(Lattice&& other) noexcept{
        std::swap(lattice, other.lattice);
        std::swap(sites, other.sites);
        std::swap(size, other.size);
        return *this;
    }

    // OPERATIOSN BETWEEN LATTICES
    // Sum
    inline Lattice& operator+=( const Lattice& other ) noexcept{
        for(int i = 0; i < sites; i++)
            lattice[i] += other.lattice[i];
        return *this;
    }
    inline Lattice& operator+=( Lattice&& other ) noexcept{
        for(int i = 0; i < sites; i++)
            lattice[i] += other.lattice[i];
        return *this;
    }

    friend Lattice operator+(Lattice lhs, const Lattice& rhs) noexcept{
        lhs += rhs;
        return lhs;
    }

    friend Lattice operator+(Lattice lhs, Lattice&& rhs) noexcept{
        lhs += rhs;
        return lhs;
    }

    // Subtraction
    inline Lattice& operator-=( const Lattice& other ) noexcept{
        for(int i = 0; i < sites; i++)
            lattice[i] -= other.lattice[i];
        return *this;
    }
    inline Lattice& operator-=( Lattice&& other ) noexcept{
        for(int i = 0; i < sites; i++)
            lattice[i] -= other.lattice[i];
        return *this;
    }

    friend Lattice operator-(Lattice lhs, const Lattice& rhs) noexcept{
        lhs -= rhs;
        return lhs;
    }

    friend Lattice operator-(Lattice lhs, Lattice&& rhs) noexcept{
        lhs -= rhs;
        return lhs;
    }

    // Multiplication
    inline Lattice& operator*=( const Lattice& other ) noexcept{
        for(int i = 0; i < sites; i++)
            lattice[i] *= other.lattice[i];
        return *this;
    }
    inline Lattice& operator*=( Lattice&& other ) noexcept{
        for(int i = 0; i < sites; i++)
            lattice[i] *= other.lattice[i];
        return *this;
    }

    friend Lattice operator*(Lattice lhs, const Lattice& rhs) noexcept{
        lhs *= rhs;
        return lhs;
    }

    friend Lattice operator*(Lattice lhs, Lattice&& rhs) noexcept{
        lhs *= rhs;
        return lhs;
    }


    // OPERATIONS WITH A SCALAR
    // Sum
    inline Lattice& operator+=( double scalar ) noexcept {
        for(int i = 0; i < sites; i++)
            lattice[i] += scalar;
        return *this;
    }
    friend Lattice operator+(Lattice lhs, double scalar) noexcept {
        lhs += scalar;
        return lhs;
    }

    // Subtraction
    inline Lattice& operator-=( double scalar ) noexcept {
        for(int i = 0; i < sites; i++)
            lattice[i] -= scalar;
        return *this;
    }
    friend Lattice operator-(Lattice lhs, double scalar) noexcept {
        lhs -= scalar;
        return lhs;
    }

    // Multiplication
    inline Lattice& operator*=( double scalar ) noexcept {
        for(int i = 0; i < sites; i++)
            lattice[i] *= scalar;
        return *this;
    }
    friend Lattice operator*(Lattice lhs, double scalar) noexcept {
        lhs *= scalar;
        return lhs;
    }
};

template <typename T>
inline Lattice<T> adj(const Lattice<T> &base){
    Lattice<T> result(base.size);
    for(int i = 0; i < result.sites; i++)
        result.at(i) = ~base.at(i);
    return std::move(result);
}

template <typename T>
inline Lattice<T> adj(Lattice<T>&& base){
    Lattice<T> result(base.size);
    for(int i = 0; i < result.sites; i++)
        result.at(i) = ~base.at(i);
    return std::move(result);
}

template <typename T>
inline Lattice<T> exp(const Lattice<T> &base){
    Lattice<T> result(base.size);
    for(int i = 0; i < result.sites; i++)
        result.at(i) = exp(base.at(i));
    return std::move(result);
}

template <typename T>
inline Lattice<T> exp(Lattice<T>&& base){
    Lattice<T> result(base.size);
    for(int i = 0; i < result.sites; i++)
        result.at(i) = exp(base.at(i));
    return std::move(result);
}


template <typename T>
inline T sum(Lattice<T>& base){
    T result = 0;
    for(int i = 0; i < base.sites; i++)
        result += base.at(i);
    return result;
}

template <>
inline SU3 sum(Lattice<SU3>& base){
    SU3 result;
    result.setSU3Zero();
    for(int i = 0; i < base.sites; i++)
        result += base.at(i);
    return result;
}

template <typename T>
inline T sum(Lattice<T>&& base){
    T result = 0;
    for(int i = 0; i < base.sites; i++)
        result += base.at(i);
    return result;
}

template <>
inline SU3 sum(Lattice<SU3>&& base){
    SU3 result;
    result.setSU3Zero();
    for(int i = 0; i < base.sites; i++)
        result += base.at(i);
    return std::move(result);
}

inline Lattice<double> realTrace(Lattice<SU3> &base){
    Lattice<double> result(base.size);
    for(int i = 0; i < base.sites; i++)
        result.at(i) = base.at(i).realTrace();
    return std::move(result);
}

inline Lattice<double> imagTrace(Lattice<SU3> &base){
    Lattice<double> result(base.size);
    for(int i = 0; i < base.sites; i++)
        result.at(i) = base.at(i).imagTrace();
    return std::move(result);
}

inline Lattice<double> realTrace(Lattice<SU3> &&base){
    Lattice<double> result(base.size);
    for(int i = 0; i < base.sites; i++)
        result.at(i) = base.at(i).realTrace();
    return std::move(result);
}

inline Lattice<double> imagTrace(Lattice<SU3> &&base){
    Lattice<double> result(base.size);
    for(int i = 0; i < base.sites; i++)
        result.at(i) = base.at(i).imagTrace();
    return std::move(result);
}


template <typename T>
Lattice<T>::Lattice(){

}


template <typename T>
Lattice<T>::Lattice(std::array<int, 4> sizeArray){
    allocate(sizeArray);
}

template <typename T>
void Lattice<T>::allocate(std::array<int, 4> sizeArray){
    // Allocate memory
    size = sizeArray;
    sites = size[0] * size[1] * size[2] * size[3];
    lattice.resize(sites);
}


void setLatticeImagIdentityValue(Lattice<SU3>& su3lat, const Lattice<double>& lat);
void setLatticeImagIdentityValue(Lattice<SU3>& su3lat, Lattice<double>&& lat);
void setToZero(Lattice<SU3>& su3lat);
void setToIdentity(Lattice<SU3>& su3lat);









template<typename T>
Lattice<T> shift(const Lattice<T>& lat, int shiftDir, int shiftStep){
    // Create Shifted Matrix
    Lattice<T> shifted;
    shifted.allocate(lat.size);

    // Create send and recv buffers
    std::array<int,4> bufferSize;
    Lattice<T> sendBuffer, recvBuffer;
    for(int i = 0; i < 4; i++)
        i != shiftDir ? bufferSize[i] = shifted.size[i] : bufferSize[i] = 1;
    sendBuffer.allocate(bufferSize);
    recvBuffer.allocate(bufferSize);

    // Fill send buffer
    // X
    if(shiftDir == 0 && shiftStep == 1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(0,y,z,t) = lat.at(0,y,z,t);
        }}}
    }
    else if(shiftDir == 0 && shiftStep == -1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(0,y,z,t) = lat.at(lat.size[0] - 1,y,z,t);
        }}}
    }

    // Y
    else if(shiftDir == 1 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,0,z,t) = lat.at(x,0,z,t);
        }}}
    }
    else if(shiftDir == 1 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,0,z,t) = lat.at(x,lat.size[1] - 1,z,t);
        }}}
    }

    // Z
    else if(shiftDir == 2 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,y,0,t) = lat.at(x,y,0,t);
        }}}
    }
    else if(shiftDir == 2 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,y,0,t) = lat.at(x,y,lat.size[2] - 1,t);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            sendBuffer.at(x,y,z,0) = lat.at(x,y,z,0);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            sendBuffer.at(x,y,z,0) = lat.at(x,y,z,lat.size[3] - 1);
        }}}
    }

    // Send non-blockingly the msg
    MPI_Request sendRequest, recvRequest;
    MPI_Status status;
    int bufferCount = 18*bufferSize[0]*bufferSize[1]*bufferSize[2]*bufferSize[3];
    int shiftRank;
    shiftStep == -1 ? shiftRank = 0 : shiftRank = 1;

    MPI_Isend(sendBuffer.lattice.data(), bufferCount, MPI_DOUBLE,
              Parallel::getNeighbor(shiftDir, abs(shiftRank-1)), 0,
              MPI_COMM_WORLD, &sendRequest);
    MPI_Irecv(recvBuffer.lattice.data(), bufferCount, MPI_DOUBLE,
              Parallel::getNeighbor(shiftDir,  shiftRank), 0,
              MPI_COMM_WORLD, &recvRequest);


    // Shift the internal points
    // X
    if(shiftDir == 0 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]-1; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x+1,y,z,t);
        }}}}
    }
    else if(shiftDir == 0 && shiftStep == -1){
        for(int x = 1; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x-1,y,z,t);
        }}}}
    }
    // Y
    else if(shiftDir == 1 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]-1; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y+1,z,t);
        }}}}
    }
    else if(shiftDir == 1 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 1; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y-1,z,t);
        }}}}
    }
    // Z
    else if(shiftDir == 2 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]-1; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z+1,t);
        }}}}
    }
    else if(shiftDir == 2 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 1; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z-1,t);
        }}}}
    }
    // T
    else if(shiftDir == 3 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]-1; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z,t+1);
        }}}}
    }
    else if(shiftDir == 3 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 1; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z,t-1);
        }}}}
    }



    // Wait for the message to arrive
    MPI_Wait(&sendRequest, &status);
    MPI_Wait(&recvRequest, &status);


    // Complete the result with the received data
    // X
    if(shiftDir == 0 && shiftStep == 1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(lat.size[0] - 1,y,z,t) = recvBuffer.at(0,y,z,t);
        }}}
    }
    else if(shiftDir == 0 && shiftStep == -1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(0,y,z,t) = recvBuffer.at(0,y,z,t);
        }}}
    }

    // Y
    else if(shiftDir == 1 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,lat.size[1] - 1,z,t) = recvBuffer.at(x,0,z,t);
        }}}
    }
    else if(shiftDir == 1 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,0,z,t) = recvBuffer.at(x,0,z,t);
        }}}
    }

    // Z
    else if(shiftDir == 2 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,lat.size[2] - 1,t) = recvBuffer.at(x,y,0,t);
        }}}
    }
    else if(shiftDir == 2 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,0,t) = recvBuffer.at(x,y,0,t);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            shifted.at(x,y,z,lat.size[3] - 1) = recvBuffer.at(x,y,z,0);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            shifted.at(x,y,z,0) = recvBuffer.at(x,y,z,0);
        }}}
    }

    // Return
    return shifted;
}


template<typename T>
Lattice<T> shift(Lattice<T>&& lat, int shiftDir, int shiftStep){
    // Create Shifted Matrix
    Lattice<T> shifted;
    shifted.allocate(lat.size);

    // Create send and recv buffers
    std::array<int,4> bufferSize;
    Lattice<T> sendBuffer, recvBuffer;
    for(int i = 0; i < 4; i++)
        i != shiftDir ? bufferSize[i] = shifted.size[i] : bufferSize[i] = 1;
    sendBuffer.allocate(bufferSize);
    recvBuffer.allocate(bufferSize);

    // Fill send buffer
    // X
    if(shiftDir == 0 && shiftStep == 1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(0,y,z,t) = lat.at(0,y,z,t);
        }}}
    }
    else if(shiftDir == 0 && shiftStep == -1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(0,y,z,t) = lat.at(lat.size[0] - 1,y,z,t);
        }}}
    }

    // Y
    else if(shiftDir == 1 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,0,z,t) = lat.at(x,0,z,t);
        }}}
    }
    else if(shiftDir == 1 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,0,z,t) = lat.at(x,lat.size[1] - 1,z,t);
        }}}
    }

    // Z
    else if(shiftDir == 2 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,y,0,t) = lat.at(x,y,0,t);
        }}}
    }
    else if(shiftDir == 2 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            sendBuffer.at(x,y,0,t) = lat.at(x,y,lat.size[2] - 1,t);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            sendBuffer.at(x,y,z,0) = lat.at(x,y,z,0);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            sendBuffer.at(x,y,z,0) = lat.at(x,y,z,lat.size[3] - 1);
        }}}
    }

    // Send non-blockingly the msg
    MPI_Request sendRequest, recvRequest;
    MPI_Status status;
    int bufferCount = 18*bufferSize[0]*bufferSize[1]*bufferSize[2]*bufferSize[3];
    int shiftRank;
    shiftStep == -1 ? shiftRank = 0 : shiftRank = 1;

    MPI_Isend(sendBuffer.lattice.data(), bufferCount, MPI_DOUBLE,
              Parallel::getNeighbor(shiftDir, abs(shiftRank-1)), 0,
              MPI_COMM_WORLD, &sendRequest);
    MPI_Irecv(recvBuffer.lattice.data(), bufferCount, MPI_DOUBLE,
              Parallel::getNeighbor(shiftDir,  shiftRank), 0,
              MPI_COMM_WORLD, &recvRequest);


    // Shift the internal points
    // X
    if(shiftDir == 0 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]-1; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x+1,y,z,t);
        }}}}
    }
    else if(shiftDir == 0 && shiftStep == -1){
        for(int x = 1; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x-1,y,z,t);
        }}}}
    }
    // Y
    else if(shiftDir == 1 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]-1; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y+1,z,t);
        }}}}
    }
    else if(shiftDir == 1 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 1; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y-1,z,t);
        }}}}
    }
    // Z
    else if(shiftDir == 2 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]-1; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z+1,t);
        }}}}
    }
    else if(shiftDir == 2 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 1; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z-1,t);
        }}}}
    }
    // T
    else if(shiftDir == 3 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]-1; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z,t+1);
        }}}}
    }
    else if(shiftDir == 3 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 1; t < lat.size[3]; t++){
            shifted.at(x,y,z,t) = lat.at(x,y,z,t-1);
        }}}}
    }



    // Wait for the message to arrive
    MPI_Wait(&sendRequest, &status);
    MPI_Wait(&recvRequest, &status);


    // Complete the result with the received data
    // X
    if(shiftDir == 0 && shiftStep == 1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(lat.size[0] - 1,y,z,t) = recvBuffer.at(0,y,z,t);
        }}}
    }
    else if(shiftDir == 0 && shiftStep == -1){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(0,y,z,t) = recvBuffer.at(0,y,z,t);
        }}}
    }

    // Y
    else if(shiftDir == 1 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,lat.size[1] - 1,z,t) = recvBuffer.at(x,0,z,t);
        }}}
    }
    else if(shiftDir == 1 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int z = 0; z < lat.size[2]; z++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,0,z,t) = recvBuffer.at(x,0,z,t);
        }}}
    }

    // Z
    else if(shiftDir == 2 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,lat.size[2] - 1,t) = recvBuffer.at(x,y,0,t);
        }}}
    }
    else if(shiftDir == 2 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int t = 0; t < lat.size[3]; t++){
            shifted.at(x,y,0,t) = recvBuffer.at(x,y,0,t);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == 1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            shifted.at(x,y,z,lat.size[3] - 1) = recvBuffer.at(x,y,z,0);
        }}}
    }
    else if(shiftDir == 3 && shiftStep == -1){
        for(int x = 0; x < lat.size[0]; x++){
        for(int y = 0; y < lat.size[1]; y++){
        for(int z = 0; z < lat.size[2]; z++){
            shifted.at(x,y,z,0) = recvBuffer.at(x,y,z,0);
        }}}
    }

    // Return
    return shifted;
}
