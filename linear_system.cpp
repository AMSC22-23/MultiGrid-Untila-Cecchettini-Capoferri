#include "linear_system.hpp"

using namespace AMG;

using namespace std;

template<typename T>
const T PoissonMatrix<T>::coeffRef(const size_t i, const size_t j){
    //these are the entries of the matrix relative to the boundary conditions; we want them to not be changed
    //by the solvers, so they will be the only entries in the whole row
    if(m_domain.isOnBoundary(m_domain.mask(i)))
        return ((j == i) ? 1. : 0.);

    if(j == i)
        return 4. * m_const_alfa / k; 
    
    //using k,l as indices on the grid
    auto [k_i, l_i] = m_domain.meshIdx(m_domain.mask(i));
    auto [k_j, l_j] = m_domain.meshIdx(m_domain.mask(j));

    size_t step = m_domain.getStep();

    // We had to set a treshold higher than h^2 due to the inexact arithmetics
    if((abs(k_i - k_j) == step) || (abs(l_i - l_j) == step))
        return - m_const_alfa / k;
    else
        return 0.;
    
}
