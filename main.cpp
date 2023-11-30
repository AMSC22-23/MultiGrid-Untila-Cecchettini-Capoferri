#include "domain.hpp"
#include "system.hpp"
#include <vector>
#include <iostream>
#include <tuple>
#include <fstream>
#include <cmath>

using namespace AMG;

double f(const double x, const double y){
    return (-5.0)*exp(x)*exp(-2.0*y);
}

double g(const double x, const double y){
    //double r = sqrt(x*x + y*y);
    return exp(x)*exp(-2.0*y);
}

template<class Vector>
void saveVectorOnFile(Vector f, std::string fileName){
    
    std::ofstream file;
    file.open(fileName);
    file<<f.size()<<std::endl;

    for(size_t i = 0; i < f.size(); i++){
        file<<f[i]<<std::endl;
    }
    file.close();
}


int main(int argc, char** argv){

int maxit= 100;
    size_t size = std::atoi(argv[1]);
  //size_t size =7;
    AMG::SquareDomain dominio(size,1.0);
    AMG::PoissonMatrix<double> A(dominio);

    AMG::DataVector<double> fVec(dominio,f,g);
    
    //saveMatrixOnFile<AMG::PoissonMatrix<double>>(A,"out.mtx");
    //saveVectorOnFile<AMG::DataVector<double>>(fVec,"f.mtx");

    std::vector<double> x(fVec.size(),0.);      //our initial guess
   
   /*
   for(int i=0;i<fVec.size();i++){
    std::cout<<x[i]<<std::endl<<std::endl;}*/

    AMG::GS<double, AMG::DataVector<double>> GAUSS(maxit, A, fVec);
    GAUSS.iteration_method(x);

    
    saveVectorOnFile<std::vector<double>>(x,"x.mtx");
    /*
    for(int i=0;i<fVec.size();i++){
     std::cout<<x[i]<<std::endl<<std::endl;}
    */

    
    return 0;
}








