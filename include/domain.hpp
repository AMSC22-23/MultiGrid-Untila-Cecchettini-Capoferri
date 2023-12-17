#ifndef DOMAIN_H
#define DOMAIN_H


#include "allIncludes.hpp"

namespace MultiGrid{

class Domain{
    public:
        virtual std::tuple<double, double> coord(const size_t i, const size_t j) const = 0;

        virtual std::tuple<size_t, size_t> meshIdx(size_t l) const = 0;
        
		//We need an access operator to the abstract domain nodes
        virtual std::tuple<double,double> operator[](const size_t i) const = 0;

		//We need a function to check if a node is on the boundary or not
        virtual bool isOnBoundary(const size_t l) const = 0;

		//to optimize the computation for each iteration of itterative solver we could need to know the non zero entries of a row
        virtual std::vector<size_t> &inRowConnections(const size_t l) = 0;

        virtual size_t mask(const size_t l) const = 0;

        virtual size_t getWidth() const = 0;

		//some useful methods
        virtual size_t numBoundaryNodes() const = 0;

        virtual size_t numConnections() const = 0;

        virtual size_t N() const = 0;

        virtual double h() const = 0;

        virtual size_t getStep() const = 0;
};


class SquareDomain: public Domain{
    private:
        size_t m_size;
        size_t step;
        size_t m_level;
        size_t width;
        double m_length;
        double m_h;
        std::vector<size_t> m_vec;

    public:

        //We need this first constructor to create a domain at a given level
        SquareDomain(const size_t size, const double length, const size_t level);

        //This second constructor creates a subdomain of the given domain, but goes down by a level
        SquareDomain(const SquareDomain &dom):SquareDomain(dom.m_size, dom.m_length, dom.m_level + 1){}

        inline std::tuple<size_t, size_t> meshIdx(size_t l) const{
            return {l / m_size, l % m_size};
        }

		//we could need an operator to get the coordinates of a specific node, for example if we have to map a function on the nodes
        inline std::tuple<double,double> coord(const size_t i, const size_t j) const {return {j* m_h, m_length - i * m_h};}
        
        std::tuple<double,double> operator[](const size_t l) const override;

        bool isOnBoundary(const size_t l) const override;

        std::vector<size_t> &inRowConnections(const size_t l) override;

        inline size_t mask(const size_t l) const{
            return step * (l / width) * m_size + step * (l % width);
        }

        inline size_t getWidth() const{return width;}

        inline size_t numBoundaryNodes() const{return width * 4 - 4;}

        inline size_t numConnections() const{return (4 * (width * width - numBoundaryNodes()));}

        inline size_t N() const{return width * width;}

        inline double h() const{return m_h * step;}

        inline size_t getStep() const{return step;}

        ~SquareDomain() = default;

};

}


#endif