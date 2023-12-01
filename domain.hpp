#include <vector>
#include <iostream>

#include <cstdlib>
//#include <Eigen/SparseCore>
//#include <unsupported/Eigen/SparseExtra>

#ifndef CLASSES_H
#define CLASSES_H

#include <tuple>
#include <vector>
#include <functional>

namespace AMG{

/*
For generalization we decided to create an abstract class Domain containing the declaration of all the public methods
we're gonna use.
*/
class Domain{
    public:
		//We need an access operator to the abstract domain nodes
        virtual std::tuple<double,double> operator[](const size_t i) const = 0;
        //getter
        virtual const size_t getWidth() const = 0;

		//We need a function to check if a node is on the boundary or not
        virtual const bool isOnBoundary(const size_t l) const = 0;

		//to optimize the computation for each iteration of itterative solver we could need to know the non zero entries of a row
        virtual const std::vector<size_t> inRowConnections(const size_t l) const = 0;

		//some useful methods
        virtual const size_t numBoundaryNodes() const = 0;
        virtual const size_t numConnections() const = 0;
        virtual const size_t N() const = 0;
        virtual const double h() const = 0;
};

/*
The goal of our project is to do some tests on a simple square domain, but the main idea was to generalize the abstract class in such
a way that we can easily make a more generic domain divided in triangles
*/
class SquareDomain: public Domain{
    public:
        SquareDomain(const size_t size, const double length):m_size(size), m_length(length){
            m_h = m_length / (m_size - 1);
        }

		//we could need an operator to get the coordinates of a specific node, for example if we have to map a function on the nodes
        std::tuple<double, double> coord(const size_t i, const size_t j) const {return {j* m_h, m_length - i * m_h};}

        std::tuple<double,double> operator[](const size_t l) const override{
            size_t i = l / m_size;
            size_t j = l % m_size;
            return coord(i,j); 
        }
       
        //getter
         const size_t getWidth() const override{
                return width;
         }


        const bool isOnBoundary(const size_t l) const override{
            size_t i = l / m_size;
            size_t j = l % m_size;
            return (((i == 0) || (j == 0) || (i == (m_size-1)) || (j == (m_size-1))) ? true : false);
        }

        const std::vector<size_t> inRowConnections(const size_t l) const override{
            std::vector<size_t> temp;
			if(isOnBoundary(l)){
				temp.push_back(l);
			}else{
				temp.push_back(l - m_size);
				temp.push_back(l - 1);
				temp.push_back(l);
				temp.push_back(l + 1);
				temp.push_back(l + m_size);
			}
			return temp;
        }

        const size_t numBoundaryNodes() const override{return m_size * 4 - 4;}
        const size_t numConnections() const override{
            return (4 * (m_size * m_size - numBoundaryNodes()));
        }
        const size_t N() const override{return m_size*m_size;}
        const double h() const override{return m_h;}

        ~SquareDomain(){}

    private:
        size_t m_size;
        double m_length;
        double m_h;
};

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

        
		
		const std::vector<size_t> nonZerosInRow(int row) const {
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


// un risolutore temporaneo

void jacobiIteration(PoissonMatrix<double> &A, DataVector<double> &f, std::vector<double> &x){
    for(size_t i = 0; i < f.size(); i++){
        double sum = 0;
        for(const auto &id : A.nonZerosInRow(i)){
            if(id != i){
                sum += A.coeffRef(i,id) * x[id];
            }
        }
        x[i] = (f[i] - sum) / A.coeffRef(i,i);
    }
}


void Interpolation(std::vector<double> &sol,Domain &domain_sup, Domain &domain_inf, DataVector<double> &f){
    for(size_t i = 0; i<domain_inf.N()-domain_inf.getWidth();i++){
        size_t index1 = domain_inf.mask(i);
        size_t index2 = domain_inf.mask(i + domain_inf.getWidth());
        size_t index3 = (index1 + index2)/2;
        //std::cout<<index3<<std::endl;
        if(! domain_sup.isOnBoundary(index3)){
            sol[index3]= 0.5*(sol[index1]+ sol[index2]);
           }else{
                sol[index3]=f[index3];
        }
    }
    size_t width = domain_sup.getWidth();
    for(size_t i=0; i<(domain_sup.N()/width); i++){
        for(size_t j=1; j<width-2;j+=2){
            sol[domain_sup.mask(j)]= 0.5*(domain_sup.mask(j-1)+domain_sup(j+1))
        }
    }
<<<<<<< Updated upstream:domain.hpp
=======

     size_t width = domain_sup.getWidth();
    for(size_t i=0; i<(domain_sup.N()/width); i++){
        for(size_t j=i*width+1; j<(i+1)*width-1;j+=2){
            
            if(! domain_sup.isOnBoundary(domain_sup.mask(j))){
            sol[domain_sup.mask(j)]= 0.5 * (sol[domain_sup.mask(j-1)]+sol[domain_sup.mask(j+1)]);
            std::cout<<j<<std::endl;
            }
            else{
                sol[domain_sup.mask(j)]=f[domain_sup.mask(j)];
                continue;
            }
        }

    }
}
>>>>>>> Stashed changes:classes/classes.hpp

}
}

#endif


