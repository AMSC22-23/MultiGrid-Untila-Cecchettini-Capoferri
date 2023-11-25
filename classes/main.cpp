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
        std::vector<size_t> row = A.nonZerosInRow(i);
        for(const auto& j : row){
            file<<i<<" "<<j<<" "<<A.coeffRef(i,j)<<std::endl;
        }
    }
    file.close();
}


int main(int argc, char** argv){
    size_t size = std::atoi(argv[1]);
    AMG::SquareDomain dominio(size,10.);
    AMG::PoissonMatrix<double> A(dominio);
    
    saveMatrixOnFile<AMG::PoissonMatrix<double>>(A,"out.mtx");
    
    return 0;
}