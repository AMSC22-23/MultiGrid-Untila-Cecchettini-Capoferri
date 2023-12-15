#ifndef LINEAR_H
#define LINEAR_H

#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <numeric>
#include <memory>

// TODO: spiegazione cosa fa
#include "multigrid.hpp"

// TODO: spiegazione cosa fa

#include "solver.hpp"
// TODO: spiegazione cosa fa
#include "utilities.hpp"
// TODO: spiegazione cosa fa
#include "linear_system.hpp"
// TODO: spiegazione cosa fa
#include "domain.hpp"

namespace MultiGrid {
// TODO: comment 
template<typename T>
class PoissonMatrix{
    private:
        Domain &m_domain;
        size_t m_size;
        T m_const_alfa;
        double k;

    public:
        PoissonMatrix(Domain &domain, const T const_alfa):m_domain(domain), m_size(domain.N()), m_const_alfa(const_alfa),
        k(domain.h() * domain.h()){}

        const T coeffRef(const size_t i, const size_t j);
		
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
            return m_domain.isOnBoundary(l);
        }

        inline const size_t rows(){return m_size;}
        inline const size_t cols(){return m_size;}

        ~PoissonMatrix() = default;
};


// TODO: comment
template<typename T>
class DataVector{
    private:
        Domain &m_domain;
        std::function<T(double,double)> m_f;
        std::function<T(double,double)> m_g;
        std::vector<T> m_vec;

    public:
        inline DataVector(Domain &domain, const std::function<T(double,double)> &f, const std::function<T(double,double)> &g) : m_domain(domain), m_f(f), m_g(g){
            for(size_t i = 0; i < m_domain.N(); i++){
                auto [x, y] = m_domain[i];
                T val = ((m_domain.isOnBoundary(i)) ? m_g : m_f)(x,y);
                
                m_vec.push_back(val);
            }
        }

        inline const T &operator[](const size_t i){ return m_vec[i]; }

        inline const size_t size(){ return m_vec.size(); }

        ~DataVector() = default;

};

};

#endif