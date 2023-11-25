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
void iteration_method(std::vector<T>& sol) override{
    //implement GS method
    std::vector<double> x_new;
    for(size_t i = 0; i<f.size(); i++){
        double sum1=0,sum2=0;
        for(int id : A.nonZerosInRow(i)){
            if(id < i)
                sum1+= A.coeffRef(i,id)*sol[id];
            else if(id > i)
                sum2+= A.coeffRef(i,id)*sol[id];
            
        }
        sol[i]= (f[i] - sum1 -sum2) / A.coeffRef(i,i);
    }
}
}//chiedi gestione var..

//jacobi method 





