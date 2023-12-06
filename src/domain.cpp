#include "domain.hpp"

using namespace AMG;


std::vector<size_t> SquareDomain::inRowConnections(const size_t l) const{
    auto equivalent_l = mask(l);
    std::vector<size_t> temp;
    if(isOnBoundary(equivalent_l)){
        temp.push_back(l);
    }else{
        temp.push_back(l - width);
        temp.push_back(l - 1);
        temp.push_back(l);
        temp.push_back(l + 1);
        temp.push_back(l + width);
    }
        return temp;
}



