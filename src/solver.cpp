#include "solver.hpp"

using namespace MultiGrid;


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


template<class Vector>
void Solver::Solve (std::vector<double> &x_k){
    size_t counter = m_maxit;
    x_k=x_k*m_res;
    while(m_res.Norm() > m_tol){
        if(counter>0){
            for(int i = 0; i < m_step; i++){
                x_k = x_k * m_it;
                counter -= 1;
            }
            x_k = x_k * m_res;
        }
        else{
            flag = 1;
            return;
        }         
    }
    flag = 0;
    return;
};




template<class Vector>

        void Residual::refresh_normalization_constant(){
            double k = 0;
            for(size_t i = 0; i < m_A.rows(); i++){
                double val = b[m_A.mask(i)];
                k += val * val;
            }
            norm_of_b = k;
        }

        

        void Residual::apply_iteration_to_vec(std::vector<double> &sol){
            norm = 0.;
            if(saveVector){
                for(size_t i = 0; i < m_A.rows(); i++){
                    size_t index = m_A.mask(i);
                    double sum = 0;
                    for(const auto &j : m_A.nonZerosInRow(i)){
                        sum += m_A.coeffRef(i, j) * sol[m_A.mask(j)];
                    }
                    double r = b[index] - sum;
                    (m_res->at(index)) = r;
                    norm += r * r;
                }
            }else{
                for(size_t i = 0; i < m_A.rows(); i++){
                    double sum = 0;
                    for(const auto &j : m_A.nonZerosInRow(i)){
                        sum += m_A.coeffRef(i, j) * sol[m_A.mask(j)];
                    }
                    double r = b[m_A.mask(i)] - sum;
                    norm += r * r;
                }
            }
        };
        
        


