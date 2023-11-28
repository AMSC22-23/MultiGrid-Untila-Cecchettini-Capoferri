#ifndef CLASSES_H
#define CLASSES_H

#include <tuple>
#include <vector>
#include <functional>



namespace Domain{

    /*
    For generalization we decided to create an abstract class Domain containing the declaration of all the public methods
    we're gonna use.
    */

    class General_Domain{
        public:
            //We need an access operator to the abstract domain nodes
            virtual std::tuple<double,double> operator[](const size_t i) const = 0;

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

    class SquareDomain: public General_Domain{
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

}



