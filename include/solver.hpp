#include "main.hpp"

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
