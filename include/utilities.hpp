#ifndef UTILITIES_H
#define UTILITIES_H

#include "main.hpp"


namespace MultiGrid
{
    // TODO: Spiegazione funzionalità metodo
    void Initialization_for_N(int argc, char** argv, unsigned int &N){}

    // TODO: Spiegazione funzionalità metodo
    template<class SpMat>
    void saveMatrixOnFile(SpMat A, std::string fileName){}

    // TODO: Spiegazione funzionalità metodo
    template<class Vector>
    void saveVectorOnFile(Vector f, std::string fileName){}

    // TODO: Spiegazione funzionalità metodo
    std::vector<double> formatVector(std::vector<double> &in, AMG::Domain &domain){}


    // TODO: Spiegazione funzionalità metodo
    inline double f(const double x, const double y){
        return -5.0 * exp(x) * exp(-2.0 * y);
    }
} // namespace UTL

#endif
