#pragma once
#include <vector>
#include <array>
#include "Math/latticemath.h"

template <typename T, int N>
class Field{
public:
    // Constructors
    Field(std::array<int, 4> size);
    inline Lattice<T>& operator[](int mu) { return m_lats[mu]; }

    int m_dimensions;
    std::array<int, 4> m_size;
private:
    std::array<Lattice<T>, N> m_lats;
};

template <typename T, int N>
Field<T,N>::Field(std::array<int, 4> size){
    // Allocate memory
    m_size = size;
    m_dimensions = N;
    for(int i = 0; i < N; i++)
        m_lats[i].allocate(m_size);
}

typedef Field<SU3,4> GluonField;
