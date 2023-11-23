#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include "classes.hpp"

using namespace std;

//some useful alias
using SpMat = Eigen::SparseMatrix<double>;
using Vec = Eigen::VectorXd;

template<class T>
void fillMatrix(T &mat){
    mat.coeffRef(0,0) = 1;

    for(auto i = 1; i < mat.cols(); i++){
        mat.coeffRef(i,i) = 4;
        mat.coeffRef(i-1,i) = -1;
        mat.coeffRef(i,i-1) = -1;
    }
}

int main(int argc, char** argv){

    //Let's create the test matrix and the test vectors
    constexpr size_t N = 10;
    SpMat A(N,N);
    A.reserve(3*N-2);
    fillMatrix(A);
    
    Vec xe = Vec::Ones(A.cols());
    Vec b = A*xe;


    
    return 0;
}