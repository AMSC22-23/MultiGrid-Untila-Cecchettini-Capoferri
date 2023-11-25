#include <iostream>
#include <tuple>
#include <fstream>
#include "classes.hpp"



int main(int argc, char** argv){
    size_t size = 7;
    AMG::SquareDomain dominio(size,10.);
    AMG::PoissonMatrix A(dominio);

    
    int count = 0;

    std::ofstream file;
    file.open("out.mtx");
    file<<A.rows()<<" "<<A.cols()<<std::endl;

    for(size_t i = 0; i < A.rows(); i++){
        for(size_t j = 0; j < A.cols(); j++){
            auto elem = A.coeffRef(i,j);
            if(elem != 0){
                file<<i<<" "<<j<<" "<<elem<<std::endl;
                count += 1;
            }
        }
    }

    std::cout<<"Non zero entries = "<<count<<std::endl;
    
    return 0;
}