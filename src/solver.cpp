#include "solver.hpp"


template<class Vector>
void Gauss_Siedel_iteration::apply_iteration_to_vec(std::vector<double> &sol) const override{
    for(size_t i = 0; i < m_A.rows(); i++){
        size_t index = m_A.mask(i);
        double sum = 0;
        for(const auto &id : m_A.nonZerosInRow(i)){
            if(id != i){
                sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
            }
        }
        sol[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
    }
};


template<class Vector>
    
void Jacobi_iteration::apply_iteration_to_vec(std::vector<double> &sol) const override{
    std::vector<double> x_new(sol.size());
    for(size_t i = 0; i < m_A.rows(); i++){
        size_t index = m_A.mask(i);
        double sum = 0;
        for(const auto &id : m_A.nonZerosInRow(i)){
            if(id != i){
                sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
            }
        }
        x_new[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
    }
    sol=x_new;
};

