 #include <iostream>
#include "domain.hpp"
#include <vector>
#include <cmath>
//#include <functions>


namespace AMG{

template<typename T,class T_VECT>
class Iteration{
    public:
    System(AMG::PoissonMatrix<T> &Mat,T_VECT &F) : mat(Mat), f(F){};
    virtual   void /*error*/ apply(std::vector<T>& sol) const=0;
    virtual   void /*error*/ apply_res(std::vector<T>& sol) const =0;
    AMG::PoissonMatrix<T> get_matrix_A (){return mat;}; 
    ~System(){};

    protected:
        T_VECT &f;
        AMG::PoissonMatrix<T> mat;
};


//gauss_saidel_method
template <typename T, class T_VECT>
class GS : public Iteration <T, T_VECT> {

public:
GS(AMG::PoissonMatrix<T>& Mat,T_VECT &F) : System<T, T_VECT>(Mat,F) {};
void apply(std::vector<T>& sol,) const override{
    
    for(size_t i = 0; i<this->f.size(); i++){
        double sum=0;
        for(auto &id : this->mat.nonZerosInRow(i)){
            if(id != i)
                sum+= this->mat.coeffRef(i,id)*sol[id];         
        }
        X_new[i]= (this->f[i] - sum1 -sum2) / this->mat.coeffRef(i,i);
        norm += std::pow(this->f[i]-)
    }
     sol=X_new;   
}

double apply_res (std::vector<T>& sol, std::vector<T>& res) override{
    double norm=0, sumres=0;//norma del residuo
    for(size_t i = 0; i<this->f.size(); i++){
        double sum=0;
        for(auto &id : this->mat.nonZerosInRow(i)){
                if(id != i){
                sum= this->mat.coeffRef(i,id)*sol[id];}

                sumres= this->mat.coeffRef(i,id)*sol[id];
        }
        sol[i]= (this->f[i]-sum) / this->mat.coeffRef(i,i);
        res[i]= f[i]-sumres;
    }
     sol=X_new;
}

private:
int iteration;

};


//jacobi method 
template <typename T, class T_VECT>
class Jacobi : public Iteration<T, T_VECT>{

public:
Jacobi(AMG::PoissonMatrix<T>& MMM,T_VECT &F) : System<T, T_VECT>(MMM,F), iteration(It) {};
void apply(std::vector<T>& sol) const override{
    
    std::vector<T> X_new(this->f.size(),0);
    
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
}

void apply_res (std::vector<T>& sol, std::vector<T>)& res override{    
 
 double norm=0, temp=0;
 std::vector<T> X_new(this->f.size(),0);
    
    for(size_t i = 0; i < this->f.size(); i++){
        double sum = 0;
        for(const auto &id : this->mat.nonZerosInRow(i)){
            if(id != i){
                sum += this->mat.coeffRef(i,id) * sol[id];
            }
        }
        double temp=(this->f[i] - sum)
        X_new[i] = temp / this->mat.coeffRef(i,i);
        res[i] = temp;
        norm +=  
    }
    sol = X_new;   
private:
int iteration;
};


template <typename T, class T_VECT>
class Multigrid : public Iteration<T, T_VECT>{
//inizializzi a zero l'errore


}


}


}