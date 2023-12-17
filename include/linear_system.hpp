#ifndef LS_H
#define LS_H

<<<<<<< HEAD

#include "allIncludes.hpp"


namespace MultiGrid{


=======
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <numeric>
#include <memory>
#include <math.h>


// TODO: spiegazione cosa fa
#include "multigrid.hpp"

// TODO: spiegazione cosa fa

#include "solver.hpp"
// TODO: spiegazione cosa fa
#include "utilities.hpp"
// TODO: spiegazione cosa fa
#include "domain.hpp"

namespace MultiGrid {
// TODO: comment 
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb
template<typename T>
class PoissonMatrix{
    public:
        //PoissonMatrix(Domain &domain, const T const_alfa);
        
        PoissonMatrix(Domain &domain, const T const_alfa):m_domain(domain), m_size(domain.N()), m_const_alfa(const_alfa),
        k(domain.h() * domain.h()){}
        
        //const T coeffRef(const size_t i, const size_t j);
        
        const T coeffRef(const size_t i, const size_t j){
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
        
		
<<<<<<< HEAD
		const std::vector<size_t> &nonZerosInRow(const size_t row){
			return m_domain.inRowConnections(row);
		}

        const size_t nonZeros(){
            return m_size + m_domain.numConnections();
        }

        const size_t mask(const size_t l){
            return m_domain.mask(l);
        }

        const size_t getWidth(){
            return m_domain.getWidth();
        }

        bool isOnBoundary(const size_t l){
=======
		inline const std::vector<size_t> &nonZerosInRow(const size_t row){
			return m_domain.inRowConnections(row);
		}

        inline const size_t nonZeros(){
            return m_size + m_domain.numConnections();
        }

        inline const size_t mask(const size_t l){
            return m_domain.mask(l);
        }

        inline size_t getWidth(){
            return m_domain.getWidth();
        }

        inline bool isOnBoundary(const size_t l){
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb
            return m_domain.isOnBoundary(l);
        }

        const size_t rows(){return m_size;}
        const size_t cols(){return m_size;}

        ~PoissonMatrix() = default;

    private:
        Domain &m_domain;
        size_t m_size;
        T m_const_alfa;
        double k;
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

        ~DataVector() = default;
<<<<<<< HEAD
    private:
        Domain &m_domain;
        std::function<T(double,double)> m_f;
        std::function<T(double,double)> m_g;
        std::vector<T> m_vec;
};




}

=======

};

};
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb

#endif