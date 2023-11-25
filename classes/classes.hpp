#ifndef CLASSES_H
#define CLASSES_H

#include <tuple>

namespace AMG{

class Domain{
    public:
        virtual std::tuple<double,double> operator[](const size_t i) const = 0;
        virtual bool isOnBoundary(size_t l) const = 0;
        virtual size_t N() const = 0;
        virtual double h() const = 0;
};


class SquareDomain: public Domain{
    public:
        SquareDomain(const size_t size, const double length):m_size(size), m_length(length){
            m_h = m_length / (m_size - 1);
        }

        std::tuple<double, double> coord(const size_t i, const size_t j) const {return {j* m_h, m_length - i * m_h};}

        std::tuple<double,double> operator[](const size_t l) const override{
            size_t i = l / m_size;
            size_t j = l % m_size;
            return coord(i,j); 
        }

        bool isOnBoundary(size_t l) const override{
            size_t i = l / m_size;
            size_t j = l % m_size;
            return (((i == 0) || (j == 0) || (i == (m_size-1)) || (j == (m_size-1))) ? true : false);
        }

        size_t N() const override{return m_size*m_size;}
        double h() const override{return m_h;}

        ~SquareDomain(){}

    private:
        size_t m_size;
        double m_length;
        double m_h;
};

class PoissonMatrix{
    public:
        PoissonMatrix(Domain &domain):m_domain(domain), m_size(domain.N()) {}

        const double coeffRef(const size_t i, const size_t j){
            if(m_domain.isOnBoundary(i)){
                return ((j == i) ? 1. : 0.);
            }else{
                if(j == i)
                    return 4.;
                
                auto [x_i, y_i] = m_domain[i];
                auto [x_j, y_j] = m_domain[j];
                auto sq_dist = (x_i - x_j) * (x_i - x_j) + (y_i - y_j) * (y_i - y_j);
                auto h = m_domain.h();

                if(sq_dist <= 1.1 * h * h)
                    return -1.;
                else
                    return 0.;
            }
        }

        const size_t rows(){return m_size;}
        const size_t cols(){return m_size;}

        ~PoissonMatrix(){}

    private:
        Domain &m_domain;
        size_t m_size;
};

}

#endif