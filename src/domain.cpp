#include "allIncludes.hpp"

//We need this first constructor to create a domain at a given level
MultiGrid::SquareDomain::SquareDomain(const size_t size, const double length, const size_t level):m_size(size),step(1) , m_level(level), width(size),
m_length(length), m_h(m_length / (m_size - 1)){
    m_vec.reserve(5);
    for(size_t i = 0; i < level; i++){
        width = (width + 1) / 2;
        step *= 2;
    }
}

std::tuple<double,double> MultiGrid::SquareDomain::operator[](const size_t l) const{
    auto [i, j] = meshIdx(mask(l));
    return coord(i,j); 
}

bool MultiGrid::SquareDomain::isOnBoundary(const size_t l) const{
    auto [i, j] = meshIdx(l);
    return (((i == 0) || (j == 0) || (i == (m_size-1)) || (j == (m_size-1))) ? true : false);
}

std::vector<size_t> &MultiGrid::SquareDomain::inRowConnections(const size_t l){
    auto equivalent_l = mask(l);
    if(isOnBoundary(equivalent_l)){
        m_vec = {l};
    }else{
        m_vec = {l - width, l - 1, l, l + 1, l + width};
    }
    return m_vec;
}






