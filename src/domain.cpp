#include "domain.hpp"

using namespace MultiGrid;

const std::vector<size_t> & SquareDomain::inRowConnections(const size_t l) const override{
    auto equivalent_l = mask(l);
    if(isOnBoundary(equivalent_l)){
        m_vec = {l};
    }else{
        m_vec = {l - width, l - 1, l, l + 1, l + width};
    }
    return m_vec;
}
