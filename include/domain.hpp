#ifndef DOMAIN_H
#define DOMAIN_H

#include "main.hpp"

using namespace std;

namespace DMN{

class Domain{
    public:
        virtual std::tuple<double, double> coord(const size_t i, const size_t j) const = 0;
        virtual std::tuple<size_t, size_t> meshIdx(size_t l) const = 0;
		//We need an access operator to the abstract domain nodes
        virtual std::tuple<double,double> operator[](const size_t i) const = 0;

		//We need a function to check if a node is on the boundary or not
        virtual const bool isOnBoundary(const size_t l) const = 0;

		//to optimize the computation for each iteration of itterative solver we could need to know the non zero entries of a row
        virtual const std::vector<size_t> &inRowConnections(const size_t l) = 0;

        virtual size_t mask(const size_t l) const = 0;
        virtual const size_t getWidth() const = 0;

		//some useful methods
        virtual const size_t numBoundaryNodes() const = 0;
        virtual const size_t numConnections() const = 0;
        virtual const size_t N() const = 0;
        virtual const double h() const = 0;
        virtual const size_t getStep() const = 0;
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
        inline SquareDomain(const size_t size, const double length, const size_t level):m_size(size),step(1) , m_level(level), width(size),
                            m_length(length), m_h(m_length / (m_size - 1)){
            m_vec.reserve(5);
            for(size_t i = 0; i < level; i++){
                width = (width + 1) / 2;
                step *= 2;
            }
        }

        inline SquareDomain(const SquareDomain &dom):SquareDomain(dom.m_size, dom.m_length, dom.m_level + 1){}


        inline std::tuple<size_t, size_t> meshIdx(size_t l) const override{  return {l / m_size, l % m_size}; }

		//we could need an operator to get the coordinates of a specific node, for example if we have to map a function on the nodes
        inline std::tuple<double, double> coord(const size_t i, const size_t j) const override {return {j* m_h, m_length - i * m_h};}

        inline std::tuple<double,double> operator[](const size_t l) const override{
            auto [i, j] = meshIdx(mask(l));
            return coord(i,j); 
        }

        inline const bool isOnBoundary(const size_t l) const override{
            auto [i, j] = meshIdx(l);
            return (((i == 0) || (j == 0) || (i == (m_size-1)) || (j == (m_size-1))) ? true : false);
        }

        const std::vector<size_t> &inRowConnections(const size_t l) override;
            
        inline size_t mask(const size_t l) const override{
            return step * (l / width) * m_size + step * (l % width);
        }

        inline size_t getWidth() const override{ return width; }

        inline size_t numBoundaryNodes() const override{return width * 4 - 4;}
        inline size_t numConnections() const override{ return (4 * (width * width - numBoundaryNodes())); }
        
        inline const size_t N() const override{return width * width;}
        inline const double h() const override{return m_h * step;}

        inline const size_t getStep() const override{return step;}

        ~SquareDomain() = default;

};

}

#endif