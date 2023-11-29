#include "domain.hpp"
#include "system.hpp"
#include <vector>
#include <iostream>
#include <tuple>
#include <fstream>
#include <cmath>

using namespace AMG;

double f(const double x, const double y){
    return 0.5;
}

double g(const double x, const double y){
    double r = sqrt(x*x + y*y);
    return sin(r);
}


int main(){

int maxit= 10;
//size_t size = std::atoi(argv[1]);
  size_t size = 30;
    AMG::SquareDomain dominio(size,10.0);
    AMG::PoissonMatrix<double> A(dominio);

    AMG::DataVector<double> fVec(dominio,f,g);
    
    //saveMatrixOnFile<AMG::PoissonMatrix<double>>(A,"out.mtx");
    //saveVectorOnFile<AMG::DataVector<double>>(fVec,"f.mtx");

    std::vector<double> x(fVec.size(),0.);      //our initial guess
   
   for(int i=0;i<fVec.size();i++){
    std::cout<<x[i]<<std::endl<<std::endl;}

    AMG::GS GAUSS(maxit, A, fVec);
    GAUSS.iteration_method(x);

    for(int i=0;i<fVec.size();i++){
     std::cout<<x[i]<<std::endl<<std::endl;}

        
    return 0;
}








