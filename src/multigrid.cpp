#include "allIncludes.hpp"

void MultiGrid::InterpolationClass::interpolate(std::vector<double> &vec){
    for(size_t i = 0; i < m_A_inf.rows() - m_A_inf.getWidth(); i++){
        size_t index1 = m_A_inf.mask(i);
        size_t index2 = m_A_inf.mask(i + m_A_inf.getWidth());
        size_t index3 = (index1 + index2) / 2;
        
        if(! m_A_sup.isOnBoundary(index3)){
            vec[index3] = 0.5 * (vec[index1] + vec[index2]);
        }else{
            vec[index3] = 0.5 * (vec[index1] + vec[index2]);
            //continue;
        }
    }
    
    size_t width = m_A_sup.getWidth();
    for(size_t i = 0; i < m_A_sup.rows() / width; i++){
        for(size_t j = i * width; j < (i + 1) * width - 1; j += 2){
            if(! m_A_sup.isOnBoundary(m_A_sup.mask(j+1)))
                vec[m_A_sup.mask(j + 1)] = 0.5 * (vec[m_A_sup.mask(j)] + vec[m_A_sup.mask(j + 2)]);
            else
                vec[m_A_sup.mask(j + 1)] = 0.5 * (vec[m_A_sup.mask(j)] + vec[m_A_sup.mask(j + 2)]);
                //continue;
        }
    }
}
