#ifndef CLASSES_H
#define CLASSES_H

#include <vector>


namespace AMG{

class SquareDomain{
    public:
        SquareDomain(const size_t &size, const double &length):m_size(size), m_length(length){
            m_h = m_length / (m_size - 1);
        }

        std::tuple<double, double> coord(const size_t &i, const size_t &j){return {j* m_h, m_length - i * m_h};}

        ~SquareDomain(){}

    private:
        size_t m_size;
        double m_length;
        double m_h;
};

}

#endif