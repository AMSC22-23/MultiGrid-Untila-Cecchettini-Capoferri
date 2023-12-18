#ifndef DOMAIN_H
#define DOMAIN_H


#include "allClasses.hpp"

namespace AMG{

class Domain{
    public:
        virtual std::tuple<double, double> coord(const size_t i, const size_t j) const = 0;

        virtual std::tuple<size_t, size_t> meshIdx(size_t l) const = 0;
        
		//We need an access operator to the abstract domain nodes
        virtual std::tuple<double,double> operator[](const size_t i) const = 0;

		//We need a function to check if a node is on the boundary or not
        virtual bool isOnBoundary(const size_t l) const = 0;

		//to optimize the computation for each iteration of itterative solver we could need to know the non zero entries of a row
        virtual const std::vector<size_t> &inRowConnections(const size_t l) = 0;

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
        SquareDomain(const size_t size, const double length, const size_t level);

        SquareDomain(const SquareDomain &dom);

        std::tuple<size_t, size_t> meshIdx(size_t l) const override;

		//we could need an operator to get the coordinates of a specific node, for example if we have to map a function on the nodes
        std::tuple<double,double> coord(const size_t i, const size_t j) const override;
        
        std::tuple<double,double> operator[](const size_t l) const override;

        bool isOnBoundary(const size_t l) const override;

        const std::vector<size_t> &inRowConnections(const size_t l) override;

        size_t mask(const size_t l) const override;

        size_t getWidth() const override;

        size_t numBoundaryNodes() const override;

        size_t numConnections() const override;

        size_t N() const override;

        double h() const override;

        size_t getStep() const override;

        ~SquareDomain() = default;


};


}


#endif