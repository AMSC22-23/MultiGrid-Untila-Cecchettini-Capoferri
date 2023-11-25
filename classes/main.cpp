#include <iostream>
#include <tuple>
#include <fstream>
#include "classes.hpp"


template<class SpMat>
void saveMatrixOnFile(SpMat A, std::string fileName){

    std::ofstream file;
    file.open(fileName);
    file<<A.rows()<<" "<<A.cols()<<" "<<A.nonZeros()<<std::endl;

    for(size_t i = 0; i < A.rows(); i++){
        for(size_t j = 0; j < A.cols(); j++){
            auto elem = A.coeffRef(i,j);
            if(elem != 0){
                file<<i<<" "<<j<<" "<<elem<<std::endl;
            }
        }
    }
    file.close();
}


int main(int argc, char** argv){
    size_t size = std::atoi(argv[1]);
    AMG::SquareDomain dominio(size,10.);
    AMG::PoissonMatrix<double> A(dominio);
    
    saveMatrixOnFile<AMG::PoissonMatrix<double>>(A,"out.mtx");

    std::cout<<"Non zero entries = "<<A.nonZeros()<<std::endl;
    
    return 0;
}