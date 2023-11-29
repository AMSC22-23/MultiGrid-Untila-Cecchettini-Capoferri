#include <iostream>
#include "domain.hpp"
#include <vector>
#include <functions>



template<typename T,class T_VECT>
class System{
    public:
    System(AMG::PoissonMatrix<T>& Mat,T_VECT &F ) : mat(Mat), f(F){};
    virtual void /*error*/ iteration_method(std::vector<T>& sol) const=0;
    AMG::PoissonMatrix get_matrix_A (){return mat;}; 
    void multigrid(){/*to define*/};
    ~System(){};

    protected:
        T_VECT f;
        AMG::PoissonMatrix mat;
};


//gauss_saidel_method

class GS : public System{

public:
GS(int It) : iteration(It) {};
void iteration_method(std::vector<T>& sol) override{
    
    std::vector<T> x_new(f.size(),0);
    for(size_t j=0; j<iteration; j++){
    for(size_t i = 0; i<f.size(); i++){
        double sum1=0,sum2=0;
        for(int id : A.nonZerosInRow(i)){
            if(id < i)
                sum1+= A.coeffRef(i,id)*X_new[id];
            else if(id > i)
                sum2+= A.coeffRef(i,id)*sol[id];
            
        }
        X_new[i]= (f[i] - sum1 -sum2) / A.coeffRef(i,i);
    }
     sol=X_new;
     X_new = X_new.assign(X_new.size(),0);
    }
};
private:
int iteration;
}//chiedi gestione var..


//jacobi method 
class Jacobi : public System{

public:
Jacobi(int It) : iteration(It) {};
void iteration_method(std::vector<T>& sol) override{
    
    std::vector<T> X_new(f.size(),0);
    for(j=0; j<iteration;j++){
    for(size_t i = 0; i < f.size(); i++){
        double sum = 0;
        for(const auto &id : A.nonZerosInRow(i)){
            if(id != i){
                sum += A.coeffRef(i,id) * x[id];
            }
        }
        X_new[i] = (f[i] - sum) / A.coeffRef(i,i);
    }
    x = X_new;
    X_new = X_new.assign(X_new.size(),0);
    }
};
private:
int iteration;
}