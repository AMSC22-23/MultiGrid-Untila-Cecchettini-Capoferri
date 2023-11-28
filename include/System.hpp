#include <tuple>
#include <vector>
#include <functional>
#include <iostream>
#include <cstdlib>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/SparseExtra>
#include "Domain.hpp"


/*
We decided for the linear system to create a class that behaves as a matrix, without allocating memory for matrix entries
We wanted to make poissonMatrix fully compatible with Eigen::SparseMatrix, but we didn't found an equivalent of
nonZerosInRow for eigen matrices
*/

template<typename T>
class PoissonMatrix{
    public:
        PoissonMatrix(Domain &domain):m_domain(domain), m_size(domain.N()) {}

        const T coeffRef(const size_t i, const size_t j){
            double k = m_domain.h() * m_domain.h();
            if(m_domain.isOnBoundary(i)){
                return ((j == i) ? 1. : 0.);
            }else{
                if(j == i)
                    return 4./k;
                
                auto [x_i, y_i] = m_domain[i];
                auto [x_j, y_j] = m_domain[j];
                auto sq_dist = (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j);
                auto h = m_domain.h();

                // We had to set a treshold higher than h^2 due to the inexact arithmetics
                if(sq_dist <= 1.1 * h * h)
                    return -1./k;
                else
                    return 0.;
            }
        }
        
        const std::vector<size_t> nonZerosInRow(const size_t row){
            return m_domain.inRowConnections(row);
        }

        const size_t nonZeros(){
            return m_size + m_domain.numConnections();
        }

        const size_t rows(){return m_size;}
        const size_t cols(){return m_size;}

        ~PoissonMatrix(){}

    private:
        Domain &m_domain;
        size_t m_size;
};

template<typename T>
class DataVector{
    public:
        DataVector(Domain &domain, const std::function<T(double,double)> &f, const std::function<T(double,double)> &g) : m_domain(domain), m_f(f), m_g(g){
            for(size_t i = 0; i < m_domain.N(); i++){
                auto [x, y] = m_domain[i];
                T val = ((m_domain.isOnBoundary(i)) ? m_g : m_f)(x,y);
                
                m_vec.push_back(val);
            }
        }
       

        const T &operator[](const size_t i){
            return m_vec[i];
        }

        const size_t size(){
            return m_vec.size();
        }

        ~DataVector(){}
    private:
        Domain &m_domain;
        std::function<T(double,double)> m_f;
        std::function<T(double,double)> m_g;
        std::vector<T> m_vec;
};

template<typename T,class T_VECT>
class System{
    public:
        System(PoissonMatrix<T>& Mat,T_VECT &F ) : mat(Mat), f(F){};
        virtual void /*error*/ iteration_method(std::vector<T>& sol) const=0;
        PoissonMatrix get_matrix_A (){return mat;}; 
        void set_matrix_A (PoissonMatrix MM){ mat = MM;}; 
        void multigrid(){/*to define*/};
        ~System(){};

    protected:
        T_VECT f;
        PoissonMatrix mat;
};