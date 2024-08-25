#ifndef SOLVERS_H
#define SOLVERS_H

#include "allIncludes.hpp"
#define TOL 1e-3

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
class Gauss_Seidel_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
    public:

        Gauss_Seidel_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {}
        
        void apply_iteration_to_vec(std::vector<double> &sol) override{
            for(size_t i = 0; i < m_A.rows(); i++){
                size_t index = m_A.mask(i); // application of multigrid
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
class BiCGSTAB : public SmootherClass<Vector> {
private:
    PoissonMatrix<double> &m_A;
    Vector &b;
    std::vector<double> r;        // Residuo
    std::vector<double> r_tilde;  // Residuo trasposto
    std::vector<double> p;        // Direzione di ricerca
    std::vector<double> v;        // Vettore intermedio
    std::vector<double> s;        // Vettore intermedio
    std::vector<double> t;        // Vettore intermedio
    double rho_old;               // Variabile per memorizzare il valore precedente di rho
    double alpha;                 // Fattore di scala per la direzione p
    double beta;                  // Fattore di scala per aggiornare p
    double omega;                 // Fattore di scala per il termine di correzione
    double tol;                   // Tolleranza per la convergenza

public:
    BiCGSTAB(PoissonMatrix<double> &A, Vector &f, double tolerance = TOL)
        : m_A(A), b(f), tol(tolerance) {
        size_t size = m_A.rows();  // Dimensione della sotto-matrice
        r.resize(size, 0.0);
        r_tilde.resize(size, 0.0);
        p.resize(size, 0.0);
        v.resize(size, 0.0);
        s.resize(size, 0.0);
        t.resize(size, 0.0);
    }

    void apply_iteration_to_vec(std::vector<double> &sol) override {
        size_t size = r.size();  // Dimensione dei vettori di lavoro
        std::cout << "Avviamento del metodo BiCGSTAB." << std::endl;

        // Calcola il residuo iniziale r = b - A * sol
        for (size_t i = 0; i < size; ++i) {
            double sum = 0.0;
            for (const auto &id : m_A.nonZerosInRow_a(i)) {
                sum += m_A.coeffRef(i, id) * sol[m_A.mask(id)];
            }
            r[i] = b[m_A.mask(i)] - sum;
        }

        // Inizializza r_tilde = r
        r_tilde = r;

        // Inizializza p a 0
        std::fill(p.begin(), p.end(), 0.0);

        // Calcola rho_0 = dot(r_tilde, r)
        double rho_0 = dot(r_tilde, r);
        rho_old = rho_0;
        omega = 1.0; // Inizializzazione di omega
        alpha = 1.0; // Inizializzazione di alpha

        // Iterazioni principali del metodo BiCGSTAB
        for (size_t k = 0; k < size; ++k) {
            // Calcola rho = dot(r_tilde, r)
            double rho = dot(r_tilde, r);

            // Calcola beta = (rho / rho_old) * (alpha / omega), evita divisione per zero
            beta = (rho_old != 0 && omega != 0) ? (rho / rho_old) * (alpha / omega) : 0.0;

            // Aggiorna p = r + beta * (p - omega * v)
            for (size_t i = 0; i < size; ++i) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            // Calcola v = A * p
            for (size_t i = 0; i < size; ++i) {
                double sum = 0.0;
                for (const auto &id : m_A.nonZerosInRow_a(i)) {
                    sum += m_A.coeffRef(i, id) * p[id];
                }
                v[i] = sum;
            }

            // Calcola alpha = rho / dot(r_tilde, v), evita divisione per zero
            double dot_r_tilde_v = dot(r_tilde, v);
            alpha = (dot_r_tilde_v != 0) ? (rho / dot_r_tilde_v) : 0.0;

            // Calcola s = r - alpha * v
            for (size_t i = 0; i < size; ++i) {
                s[i] = r[i] - alpha * v[i];
            }

            // Calcola t = A * s
            for (size_t i = 0; i < size; ++i) {
                double sum = 0.0;
                for (const auto &id : m_A.nonZerosInRow_a(i)) {
                    sum += m_A.coeffRef(i, id) * s[id];
                }
                t[i] = sum;
            }

            // Calcola omega = dot(t, s) / dot(t, t), evita divisione per zero
            double dot_t_t = dot(t, t);
            omega = (dot_t_t != 0) ? (dot(t, s) / dot_t_t) : 0.0;

            // Aggiorna la soluzione sol = sol + alpha * p + omega * s
            for (size_t i = 0; i < size; ++i) {
                size_t sol_index = m_A.mask(i); // usa mask per mappare l'indice
                sol[sol_index] += alpha * p[i] + omega * s[i];
            }

            // Aggiorna il residuo r = s - omega * t
            for (size_t i = 0; i < size; ++i) {
                r[i] = s[i] - omega * t[i];
            }

            // Calcola rho_old per la prossima iterazione
            rho_old = rho;

            // Controlla la condizione di arresto
            double residual_norm = std::sqrt(dot(r, r));
            std::cout << "Norma del residuo: " << residual_norm << std::endl;
            if (residual_norm < tol) {
                std::cout << "Convergenza raggiunta." << std::endl;
                break;
            }
        }
    }

private:
    double dot(const std::vector<double> &v1, const std::vector<double> &v2) const {
        double sum = 0.0;
        for (size_t i = 0; i < v1.size(); ++i) {
            sum += v1[i] * v2[i];
        }
        return sum;
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