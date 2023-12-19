#ifndef SOLVERS_H
#define SOLVERS_H

#include "allIncludes.hpp"

/*// versione ancora da ridefinire e testare
template<class Matrix, class Vector>
void apply_iteration_to_vec_parallel(std::vector<double> &sol, Matrix m_A, Vector b){
    std::vector<double> x_new(sol.size());

        #pragma omp parallel for
        for(size_t i = 0; i < m_A.rows(); i++){
            size_t index = m_A.mask(i);
            double sum = 0;

            for(const auto &id : m_A.nonZerosInRow(i)){
                if(id != i){
                    sum += m_A.coeffRef(i, id) * sol[m_A.mask(id)];
                }
            }
            x_new[index] = (b[index] - sum) / m_A.coeffRef(i, i);
        }
        sol = x_new;
    
}*/

namespace MultiGrid{

template<class Vector>
class SmootherClass{
    public:

        virtual void apply_iteration_to_vec(std::vector<double> &sol) = 0;

        inline friend std::vector<double>& operator*(std::vector<double> &x_k, SmootherClass &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        } // x^(k+1) = x^(k) * B

};


template<class Vector>
class Gauss_Siedel_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
    public:

        Gauss_Siedel_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {}
        
        void apply_iteration_to_vec(std::vector<double> &sol) override{
            for(size_t i = 0; i < m_A.rows(); i++){
                size_t index = m_A.mask(i);
                double sum = 0;
                if(m_A.isOnBoundary(m_A.mask(i))){
                    sum = 0;
                }else{
                    for(const auto &id : m_A.nonZerosInRow_a(i)){
                        if(id != i){
                            sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
                        }
                    }
                }
                sol[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
            }
        }
};



template<class Vector>
class Jacobi_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
        std::vector<double> temp;
    public:

        Jacobi_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {
            temp = std::vector<double>(b.size(),0.);
        }
       
        void apply_iteration_to_vec(std::vector<double> &sol) override{
            for(size_t i = 0; i < m_A.rows(); i++){
                size_t index = m_A.mask(i);
                double sum = 0;
                if(m_A.isOnBoundary(m_A.mask(i))){
                    sum = 0;
                }else{
                    for(const auto &id : m_A.nonZerosInRow_a(i)){
                        if(id != i){
                            sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
                        }
                    }
                }
                temp[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
            }
            sol.swap(temp);
        }
// one iteration of jacobi
};

template<class Vector>
class Parallel_Jacobi_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
        std::vector<double> temp;
    public:

        Parallel_Jacobi_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {
            temp = std::vector<double>(b.size(),0.);
        }
        void apply_iteration_to_vec(std::vector<double> &sol) override{
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for(size_t i = 0; i < m_A.rows(); i++){
                size_t index = m_A.mask(i);
                double sum = 0;
                if(m_A.isOnBoundary(m_A.mask(i))){
                    sum = 0;
                }else{
                    for(const auto &id : m_A.nonZerosInRow_a(i)){
                        if(id != i){
                            sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
                        }
                    }
                }
                temp[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
            }
            sol.swap(temp);
        }
};


template<class Vector>
class Residual{
    private:
        PoissonMatrix<double> &m_A;
        Vector &b;
        std::vector<double> *m_res;
        bool saveVector;
        double norm_of_b;
        double norm;
        
    public:
        Residual(PoissonMatrix<double> &A, Vector &f): m_A(A), b(f), m_res(NULL), saveVector(false), norm_of_b(0.){
            for(size_t i = 0; i < b.size(); i++){
                double val = b[i];
                norm_of_b += val * val;
            }
        }

        Residual(PoissonMatrix<double> &A, Vector &f, std::vector<double> &res): m_A(A), b(f), m_res(&res), saveVector(true), norm_of_b(0.){
            for(size_t i = 0; i < A.rows(); i++){
                double val = b[A.mask(i)];
                norm_of_b += val * val;
            }
        }

        void refresh_normalization_constant(){
            double k = 0;
            #ifdef _OPENMP
            #pragma omp parallel for reduction(+:k)
            #endif
            for(size_t i = 0; i < m_A.rows(); i++){
                double val = b[m_A.mask(i)];
                k += val * val;
            }
            norm_of_b = k;
        }     


        void apply_iteration_to_vec(std::vector<double> &sol){
            norm = 0.;
            if(saveVector){
                #ifdef _OPENMP
                #pragma omp parallel for reduction(+:norm)
                #endif
                for(size_t i = 0; i < m_A.rows(); i++){
                    size_t index = m_A.mask(i);
                    double sum = 0;
                    if(m_A.isOnBoundary(index)){
                        sum = m_A.coeffRef(i,i) * sol[index];
                    }else{
                        for(const auto &j : m_A.nonZerosInRow_a(i)){
                            sum += m_A.coeffRef(i, j) * sol[m_A.mask(j)];
                        }
                    }
                    double r = b[index] - sum;
                    (m_res->at(index)) = r;
                    norm += r * r;
                }
            }else{
                #ifdef _OPENMP
                #pragma omp parallel for reduction(+:norm)
                #endif
                for(size_t i = 0; i < m_A.rows(); i++){
                    size_t index = m_A.mask(i);
                    double sum = 0;
                    if(m_A.isOnBoundary(index)){
                        sum = m_A.coeffRef(i,i) * sol[index];
                    }else{
                        for(const auto &j : m_A.nonZerosInRow_a(i)){
                            sum += m_A.coeffRef(i, j) * sol[m_A.mask(j)];
                        }
                    }
                    double r = b[index] - sum;
                    norm += r * r;
                }
            }
        }

        
        friend std::vector<double>& operator*(std::vector<double> &x_k, Residual &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }
        

        double Norm(){
            return sqrt(norm/norm_of_b);
        }
};

template<class Vector>
class Solver{

    private:
        SmootherClass<Vector> &m_it;
        Residual<Vector> &m_res;
        size_t m_maxit;
        double m_tol;
        int flag;
        int m_step;

    public:
        Solver(SmootherClass<Vector> &it,Residual<Vector> &res, size_t maxit, double tol, int step) : m_it(it), m_res(res), m_maxit(maxit), m_tol(tol), m_step(step) {};

        void Solve (std::vector<double> &x_k){
            size_t counter = m_maxit;
            x_k=x_k*m_res;
            while(m_res.Norm() > m_tol){
                if(counter>0){
                    for(int i = 0; i < m_step; i++){
                        x_k = x_k * m_it;
                        counter -= 1;
                    }
                    x_k = x_k * m_res;
                }
                else{
                    flag = 1;
                    return;
                }         
            }
            flag = 0;
            return;
        }
        int Status(){
            return flag;
        }
        
        friend std::vector<double>& operator*(std::vector<double> &x_k, Solver &B)
        {
            B.Solve(x_k);
            return x_k;
        }

};

}




#endif