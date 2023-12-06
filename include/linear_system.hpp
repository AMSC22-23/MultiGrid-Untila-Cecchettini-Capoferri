#ifndef LINEAR_H
#define LINEAR_H

#include "domain.hpp"
#include <cmath>
#include <iostream>

namespace AMG{

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
		
		inline const std::vector<size_t> nonZerosInRow(const size_t row){ return m_domain.inRowConnections(row); }

        inline const size_t nonZeros(){ return m_size + m_domain.numConnections(); }

        inline const size_t mask(const size_t l){ return m_domain.mask(l); }

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

// TODO: comment
template<class Vector>
class Iteration{
    protected:
        Vector &b; // Ax = b
    public:
        inline Iteration(Vector &sol) : b(sol) {};

        virtual void apply_iteration_to_vec(std::vector<double> &sol) const = 0;

        inline friend std::vector<double>& operator*(std::vector<double> &x_k, const Iteration &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        } // B is iteration matrix

        // return the norm of the residual
        virtual double apply_with_residual(std::vector<double> &sol) const = 0;
};



template<class Vector>
class Gauss_Siedel_iteration : public Iteration<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
    public:

        Gauss_Siedel_iteration(PoissonMatrix<double> &A, Vector &sol) : m_A(A), Iteration<Vector>(sol) {};

       
        void apply_iteration_to_vec(std::vector<double> &sol) const override;
         // TODO implementation



// one iteration of GS

};








}

#endif