#include "linear_system.hpp"

//implementation of Gauss-Seidel method

template <typename T> 
class GS : public System{

    public:
        void iteration_method(std::vector<T>& sol) override{
        
            std::vector<T> x_sol;
            double sum1,sum2;

            for(size_t i = 0; i<f.size(); i++){
                
                sum1=0;
                sum2=0;
                
                for(int idx : A.nonZerosInRow(i)){
                    if(idx < i)
                        sum1+= A.coeffRef(i,idx)*x_sol[idx];
                    else if(idx > i)
                        sum2+= A.coeffRef(i,idx)*sol[idx];
                }
                x_sol[i]= (f[i] - sum1 -sum2) / A.coeffRef(i,i);
            }
            sol=x_sol;

        // TODO: adding error norm control and tolerance 
        }
}//chiedi gestione var..

//jacobi method 
class Jacobi : public System{

    public:

        void iteration_method(std::vector<T>& sol) override{

            std::vector<T> s_sol;
            double sum;

            for(size_t i = 0; i < f.size(); i++){
                sum = 0;
                for(const auto &id : A.nonZerosInRow(i)){
                    if(id != i){
                        sum += A.coeffRef(i,id) * sol[id];
                    }
                }
                s_sol[i] = (f[i] - sum) / A.coeffRef(i,i);
            }
            sol = s_sol;
        
        }
}


// Multigrid method implementation
template <typename T, class Iterative_method, class Matrix> 
class Multigrid : public System{

    public:
        Multigrid(Iterative_method Method, unsigned int N_iteration_for_each_step) : Method(Method), 
        N_iteration_for_each_step(N_iteration_for_each_step), mat_restricted(Method.get_matrix_A()) {};

        void Restriction_mat(Matrix & Matrix_out, unsigned int & numbers_of_restrictions){
            Numbers_of_restrictions = numbers_of_restrictions;  // saved into private var number of resize attuacted

            while(numbers_of_restrictions != 0){
                // TODO: N iteration of Iterative_method
                // TODO: restriction of course of matrix A into small course grid
                numbers_of_restrictions--;
            }

        }

        void Increasing_mat(Matrix & Matrix_out){

            while(Numbers_of_restrictions != 0){
                // TODO: increase of course of matrix A into large course grid
                // QUESTION: adding iterative method also here?
                numbers_of_restrictions--;
            }
        }

    protected:
        unsigned int N_iteration_for_each_step = 10; // 10 is by default
        Iterative_method Method;
        PoissonMatrix mat_restricted;
        unsigned int Numbers_of_restrictions = 0;
}