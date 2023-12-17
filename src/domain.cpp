#include "allIncludes.hpp"

<<<<<<< HEAD
//We need this first constructor to create a domain at a given level
MultiGrid::SquareDomain::SquareDomain(const size_t size, const double length, const size_t level):m_size(size),step(1) , m_level(level), width(size),
m_length(length), m_h(m_length / (m_size - 1)){
    m_vec.reserve(5);
    for(size_t i = 0; i < level; i++){
        width = (width + 1) / 2;
        step *= 2;
    }
}

//This second constructor creates a subdomain of the given domain, but goes down by a level
MultiGrid::SquareDomain::SquareDomain(const SquareDomain &dom):SquareDomain(dom.m_size, dom.m_length, dom.m_level + 1){}

inline std::tuple<size_t, size_t> MultiGrid::SquareDomain::meshIdx(size_t l) const{
    return {l / m_size, l % m_size};
}

inline std::tuple<double,double> MultiGrid::SquareDomain::coord(const size_t i, const size_t j) const {return {j* m_h, m_length - i * m_h};}

inline std::tuple<double,double> MultiGrid::SquareDomain::operator[](const size_t l) const{
    auto [i, j] = meshIdx(mask(l));
    return coord(i,j); 
}

inline const bool MultiGrid::SquareDomain::isOnBoundary(const size_t l) const{
    auto [i, j] = meshIdx(l);
    return (((i == 0) || (j == 0) || (i == (m_size-1)) || (j == (m_size-1))) ? true : false);
}

inline const std::vector<size_t> &MultiGrid::SquareDomain::inRowConnections(const size_t l){
=======
using namespace MultiGrid;

const std::vector<size_t> & SquareDomain::inRowConnections(const size_t l) const override{
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb
    auto equivalent_l = mask(l);
    if(isOnBoundary(equivalent_l)){
        m_vec = {l};
    }else{
        m_vec = {l - width, l - 1, l, l + 1, l + width};
    }
    return m_vec;
<<<<<<< HEAD
}

inline size_t MultiGrid::SquareDomain::mask(const size_t l) const{
    return step * (l / width) * m_size + step * (l % width);
}

inline const size_t MultiGrid::SquareDomain::getWidth() const{return width;}

inline const size_t MultiGrid::SquareDomain::numBoundaryNodes() const{return width * 4 - 4;}

inline const size_t MultiGrid::SquareDomain::numConnections() const{return (4 * (width * width - numBoundaryNodes()));}

inline const size_t MultiGrid::SquareDomain::N() const{return width * width;}

inline const double MultiGrid::SquareDomain::h() const{return m_h * step;}

inline const size_t MultiGrid::SquareDomain::getStep() const{return step;}


=======
}
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb
