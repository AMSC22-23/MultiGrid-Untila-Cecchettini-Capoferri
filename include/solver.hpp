#ifndef SOLVER_H
#define SOLVER_H

#include "main.hpp"


namespace MultiGrid {

// TODO: specific of what class does
template<class Vector>
class SmootherClass{
    public:
        
        virtual void apply_iteration_to_vec(std::vector<double> &sol) const = 0;

        inline friend std::vector<double>& operator*(std::vector<double> &x_k, const SmootherClass &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        } 
};


template<class Vector>
class Gauss_Siedel_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
    public:

        Gauss_Siedel_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {}  
        void apply_iteration_to_vec(std::vector<double> &sol) const override;
};


template<class Vector>
class Jacobi_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
    public:

        Jacobi_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {};
       
        void apply_iteration_to_vec(std::vector<double> &sol) const override;
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

        void Solve (std::vector<double> &x_k)
        inline int Status(){
            return flag;
        }
        
        friend std::vector<double>& operator*(std::vector<double> &x_k, Solver &B)
        {
            B.Solve(x_k);
            return x_k;
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
        void refresh_normalization_constant();

        void apply_iteration_to_vec(std::vector<double> &sol);
        
        friend std::vector<double>& operator*(std::vector<double> &x_k, Residual &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }
        

        inline double Norm(){
            return sqrt(norm/norm_of_b);
        }
    
};

};
#endif