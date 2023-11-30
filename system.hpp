 #include <iostream>
#include "domain.hpp"
#include <vector>
//#include <functions>


namespace AMG{

template<typename T,class T_VECT>
class System{
    public:
    System(AMG::PoissonMatrix<T> &Mat,T_VECT &F) : mat(Mat), f(F){};
    void /*error*/ iteration_method(std::vector<T>& sol);
    AMG::PoissonMatrix<T> get_matrix_A (){return mat;}; 
    ~System(){};

    protected:
        T_VECT &f;
        AMG::PoissonMatrix<T> mat;
};


//gauss_saidel_method
template <typename T, class T_VECT>
class GS : public System <T, T_VECT> {

public:
GS(int It, AMG::PoissonMatrix<T>& Mat,T_VECT &F) : System<T, T_VECT>(Mat,F), iteration(It) {};
void iteration_method(std::vector<T>& sol) {
    
    std::vector<T> X_new(this->f.size(),0);
    
    for(size_t k=0; k<iteration; k++){
    for(size_t i = 0; i<this->f.size(); i++){
        double sum1=0,sum2=0;
        for(auto &id : this->mat.nonZerosInRow(i)){
            if(id < i)
                sum1+= this->mat.coeffRef(i,id)*X_new[id];
            else if(id > i)
                sum2+= this->mat.coeffRef(i,id)*sol[id];
            
        }
        X_new[i]= (this->f[i] - sum1 -sum2) / this->mat.coeffRef(i,i);
    }
     sol=X_new;
     X_new.assign(X_new.size(),0);
    }
}
private:
int iteration;

};


//jacobi method 
template <typename T, class T_VECT>
class Jacobi : public System<T, T_VECT>{

public:
Jacobi(int It, AMG::PoissonMatrix<T>& MMM,T_VECT &F) : System<T, T_VECT>(MMM,F), iteration(It) {};
void iteration_method(std::vector<T>& sol) {
    
    std::vector<T> X_new(this->f.size(),0);
    for(int k=0; k<iteration;k++){
    for(size_t i = 0; i < this->f.size(); i++){
        double sum = 0;
        for(const auto &id : this->mat.nonZerosInRow(i)){
            if(id != i){
                sum += this->mat.coeffRef(i,id) * sol[id];
            }
        }
        X_new[i] = (this->f[i] - sum) / this->mat.coeffRef(i,i);
    }
    sol = X_new;
    X_new.assign(X_new.size(),0);
    }
}
private:
int iteration;
};
}