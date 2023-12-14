#ifndef CLASSES_H
#define CLASSES_H

#include <tuple>
#include <vector>
#include <functional>
#include <numeric>
#include <memory>
#include "utils.hpp"


namespace AMG{

/*
For generalization we decided to create an abstract class Domain containing the declaration of all the public methods
we're gonna use.
*/
class Domain{
    public:
        virtual std::tuple<double, double> coord(const size_t i, const size_t j) const = 0;
        virtual std::tuple<size_t, size_t> meshIdx(size_t l) const = 0;
		//We need an access operator to the abstract domain nodes
        virtual std::tuple<double,double> operator[](const size_t i) const = 0;

		//We need a function to check if a node is on the boundary or not
        virtual const bool isOnBoundary(const size_t l) const = 0;

		//to optimize the computation for each iteration of itterative solver we could need to know the non zero entries of a row
        virtual const std::vector<size_t> &inRowConnections(const size_t l) = 0;

        virtual size_t mask(const size_t l) const = 0;
        virtual const size_t getWidth() const = 0;

		//some useful methods
        virtual const size_t numBoundaryNodes() const = 0;
        virtual const size_t numConnections() const = 0;
        virtual const size_t N() const = 0;
        virtual const double h() const = 0;
        virtual const size_t getStep() const = 0;
};

/*
The goal of our project is to do some tests on a simple square domain, but the main idea was to generalize the abstract class in such
a way that we can easily make a more generic domain divided in triangles
*/
class SquareDomain: public Domain{
    public:
        SquareDomain(const size_t size, const double length, const size_t level):m_size(size),step(1) , m_level(level), width(size),
        m_length(length), m_h(m_length / (m_size - 1)){
            m_vec.reserve(5);
            for(size_t i = 0; i < level; i++){
                width = (width + 1) / 2;
                step *= 2;
            }
        }

        SquareDomain(const SquareDomain &dom):SquareDomain(dom.m_size, dom.m_length, dom.m_level + 1){}

        std::tuple<size_t, size_t> meshIdx(size_t l) const override{
            return {l / m_size, l % m_size};
        }

		//we could need an operator to get the coordinates of a specific node, for example if we have to map a function on the nodes
        std::tuple<double, double> coord(const size_t i, const size_t j) const override {return {j* m_h, m_length - i * m_h};}

        std::tuple<double,double> operator[](const size_t l) const override{
            auto [i, j] = meshIdx(mask(l));
            return coord(i,j); 
        }

        const bool isOnBoundary(const size_t l) const override{
            auto [i, j] = meshIdx(l);
            return (((i == 0) || (j == 0) || (i == (m_size-1)) || (j == (m_size-1))) ? true : false);
        }

        const std::vector<size_t> &inRowConnections(const size_t l) override{
            auto equivalent_l = mask(l);
			if(isOnBoundary(equivalent_l)){
				m_vec = {l};
			}else{
                /*
				temp.push_back(l - width);
				temp.push_back(l - 1);
				temp.push_back(l);
				temp.push_back(l + 1);
				temp.push_back(l + width);
                */
               m_vec = {l - width, l - 1, l, l + 1, l + width};
			}
			return m_vec;
        }

        size_t mask(const size_t l) const override{
            return step * (l / width) * m_size + step * (l % width);
        }

        const size_t getWidth() const override{
            return width;
        }

        const size_t numBoundaryNodes() const override{return width * 4 - 4;}
        const size_t numConnections() const override{
            return (4 * (width * width - numBoundaryNodes()));
        }
        const size_t N() const override{return width * width;}
        const double h() const override{return m_h * step;}

        const size_t getStep() const override{return step;}

        ~SquareDomain() = default;

    private:
        size_t m_size;
        size_t step;
        size_t m_level;
        size_t width;
        double m_length;
        double m_h;
        std::vector<size_t> m_vec;

};

/*
We decided for the linear system to create a class that behaves as a matrix, without allocating memory for matrix entries
We wanted to make poissonMatrix fully compatible with Eigen::SparseMatrix, but we didn't found an equivalent of
nonZerosInRow for eigen matrices
*/

template<typename T>
class PoissonMatrix{
    public:
        PoissonMatrix(Domain &domain, const T const_alfa):m_domain(domain), m_size(domain.N()), m_const_alfa(const_alfa),
        k(domain.h() * domain.h()){}

        const T coeffRef(const size_t i, const size_t j){
            //these are the entries of the matrix relative to the boundary conditions; we want them to not be changed
            //by the solvers, so they will be the only entries in the whole row
            if(m_domain.isOnBoundary(m_domain.mask(i)))
                return ((j == i) ? 1. : 0.);

            if(j == i)
                return 4. * m_const_alfa / k; 
            
            //using k,l as indices on the grid
            auto [k_i, l_i] = m_domain.meshIdx(m_domain.mask(i));
            auto [k_j, l_j] = m_domain.meshIdx(m_domain.mask(j));

            size_t step = m_domain.getStep();

            // We had to set a treshold higher than h^2 due to the inexact arithmetics
            if((abs(k_i - k_j) == step) || (abs(l_i - l_j) == step))
                return - m_const_alfa / k;
            else
                return 0.;
            
        }
		
		const std::vector<size_t> &nonZerosInRow(const size_t row){
			return m_domain.inRowConnections(row);
		}

        const size_t nonZeros(){
            return m_size + m_domain.numConnections();
        }

        const size_t mask(const size_t l){
            return m_domain.mask(l);
        }

        const size_t getWidth(){
            return m_domain.getWidth();
        }

        bool isOnBoundary(const size_t l){
            return m_domain.isOnBoundary(l);
        }

        const size_t rows(){return m_size;}
        const size_t cols(){return m_size;}

        ~PoissonMatrix() = default;

    private:
        Domain &m_domain;
        size_t m_size;
        T m_const_alfa;
        double k;
};

template<typename T>
class DataVector{
    public:
        DataVector(Domain &domain, const std::function<T(double,double)> &f, const std::function<T(double,double)> &g) : m_domain(domain), m_f(f), m_g(g){
            for(size_t i = 0; i < m_domain.N(); i++){
                auto [x, y] = m_domain[i];
                T val = ((m_domain.isOnBoundary(i)) ? m_g : m_f)(x,y);
                
                m_vec.push_back(val);
            }
        }
       

        const T &operator[](const size_t i){
            return m_vec[i];
        }

        const size_t size(){
            return m_vec.size();
        }

        ~DataVector() = default;
    private:
        Domain &m_domain;
        std::function<T(double,double)> m_f;
        std::function<T(double,double)> m_g;
        std::vector<T> m_vec;
};

/*
// un risolutore temporaneo

template<class Vector>
void gaussSeidelIteration(PoissonMatrix<double> &A, Vector &f, std::vector<double> &x){
    for(size_t i = 0; i < A.rows(); i++){
        size_t index = A.mask(i);
        double sum = 0;
        for(const auto &id : A.nonZerosInRow(i)){
            if(id != i){
                sum += A.coeffRef(i,id) * x[A.mask(id)];
            }
        }
        x[index] = (f[index] - sum) / A.coeffRef(i,i);
    }
}


void Interpolation(std::vector<double> &sol, Domain &domain_sup, Domain &domain_inf){
    for(size_t i = 0; i < domain_inf.N() - domain_inf.getWidth(); i++){
        size_t index1 = domain_inf.mask(i);
        size_t index2 = domain_inf.mask(i + domain_inf.getWidth());
        size_t index3 = (index1 + index2) / 2;
        
        if(! domain_sup.isOnBoundary(index3)){
            sol[index3] = 0.5 * (sol[index1] + sol[index2]);
        }else{
            //sol[index3] = f[index3];
            continue;
        }
    }
    
    size_t width = domain_sup.getWidth();
    for(size_t i = 0; i < domain_sup.N() / width; i++){
        for(size_t j = i * width; j < (i + 1) * width - 1; j += 2){
            if(! domain_sup.isOnBoundary(domain_sup.mask(j+1)))
                sol[domain_sup.mask(j + 1)] = 0.5 * (sol[domain_sup.mask(j)] + sol[domain_sup.mask(j + 2)]);
            else
                //sol[domain_sup.mask(j + 1)] = f[domain_sup.mask(j + 1)];
                continue;
        }
    }
    
}

double error(std::vector<double> &vec, std::vector<double> &sol, Domain &domain){
    double squaredError = 0;
    double squaredNorm = 0;
    
    for(size_t i = 0; i < domain.N(); i++){
        size_t index = domain.mask(i);
        double solVal = sol[index];
        double errVal = vec[index] - solVal;
        squaredError += errVal * errVal;
        squaredNorm += solVal * solVal;
    }
    return sqrt(squaredError / squaredNorm); 
}


template<class Vector>
void residual(std::vector<double> &u, Vector &f, PoissonMatrix<double> &A, std::vector<double> &res){
    for(size_t i = 0; i < A.rows(); i++){
        size_t index = A.mask(i);
        double sum = 0;
        for(const auto &j : A.nonZerosInRow(i)){
            sum += A.coeffRef(i, j) * u[A.mask(j)];
        }
        res[index] = f[index] - sum;
    }
}

template<class Vector>
double residualNorm(std::vector<double> &u, Vector &f, PoissonMatrix<double> &A){
    double res = 0.;
    for(size_t i = 0; i < A.rows(); i++){
        size_t index = A.mask(i);
        double sum = 0.;
        for(const auto &j : A.nonZerosInRow(i)){
            sum += A.coeffRef(i, j) * u[A.mask(j)];
        }
        sum = f[index] - sum;
        res += sum * sum;
    }
    return sqrt(res);
}

template<class Vector>
double norm(Vector &u, PoissonMatrix<double> &A){
    double sum = 0;
    for(size_t i = 0; i < A.rows(); i++){
        double val = u[A.mask(i)];
        sum += val * val;
    }
    return sqrt(sum);
}

*/

// TODO: comment
template<class Vector>
class SmootherClass{
    protected:
        //Vector &b; // Ax = b
    public:
        //inline Iteration(Vector &sol) : b(sol) {};

        virtual void apply_iteration_to_vec(std::vector<double> &sol) const = 0;

        inline friend std::vector<double>& operator*(std::vector<double> &x_k, const SmootherClass &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        } // B is equivalent to iteration matrix

        //virtual double Res_calc(std::vector<double> &sol) const = 0;

        // return the norm of the residual
        //virtual double apply_with_residual(std::vector<double> &sol) const = 0;
};



template<class Vector>
class Gauss_Siedel_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
    public:

        Gauss_Siedel_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {}
        //Iteration.set(      
        void apply_iteration_to_vec(std::vector<double> &sol) const override{
            for(size_t i = 0; i < m_A.rows(); i++){
                size_t index = m_A.mask(i);
                double sum = 0;
                for(const auto &id : m_A.nonZerosInRow(i)){
                    if(id != i){
                        sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
                    }
                }
                sol[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
            }
        }
         // TODO implementation

        /*
        double Res_calc(std::vector<double> &sol) const override {
            std::accumulate
        }
        */
// one iteration of GS

};

template<class Vector>
class Jacobi_iteration : public SmootherClass<Vector>{
    private:    
        PoissonMatrix<double> &m_A;
        Vector &b; // Ax = b
    public:

        Jacobi_iteration(PoissonMatrix<double> &A, Vector &f) : m_A(A), b(f) {};
       
        void apply_iteration_to_vec(std::vector<double> &sol) const override{
            std::vector<double> x_new(sol.size());
            for(size_t i = 0; i < m_A.rows(); i++){
                size_t index = m_A.mask(i);
                double sum = 0;
                for(const auto &id : m_A.nonZerosInRow(i)){
                    if(id != i){
                        sum += m_A.coeffRef(i,id) * sol[m_A.mask(id)];
                    }
                }
                x_new[index] = (this->b[index] - sum) / m_A.coeffRef(i,i);
            }
            sol=x_new;
        }
         // TODO implementation


// one iteration of GS

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
            for(size_t i = 0; i < m_A.rows(); i++){
                double val = b[m_A.mask(i)];
                k += val * val;
            }
            norm_of_b = k;
        }

        

        void apply_iteration_to_vec(std::vector<double> &sol){
            norm = 0.;
            if(saveVector){
                for(size_t i = 0; i < m_A.rows(); i++){
                    size_t index = m_A.mask(i);
                    double sum = 0;
                    for(const auto &j : m_A.nonZerosInRow(i)){
                        sum += m_A.coeffRef(i, j) * sol[m_A.mask(j)];
                    }
                    double r = b[index] - sum;
                    (m_res->at(index)) = r;
                    norm += r * r;
                }
            }else{
                for(size_t i = 0; i < m_A.rows(); i++){
                    double sum = 0;
                    for(const auto &j : m_A.nonZerosInRow(i)){
                        sum += m_A.coeffRef(i, j) * sol[m_A.mask(j)];
                    }
                    double r = b[m_A.mask(i)] - sum;
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
        void debug(){
            std::cout<<norm_of_b<<std::endl;
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

class InterpolationClass{
    private:
        PoissonMatrix<double> &m_A_inf, &m_A_sup;

    public:
        InterpolationClass(PoissonMatrix<double> &A_inf, PoissonMatrix<double> &A_sup): m_A_inf(A_inf), m_A_sup(A_sup){}

        void interpolate(std::vector<double> &vec){
            for(size_t i = 0; i < m_A_inf.rows() - m_A_inf.getWidth(); i++){
                size_t index1 = m_A_inf.mask(i);
                size_t index2 = m_A_inf.mask(i + m_A_inf.getWidth());
                size_t index3 = (index1 + index2) / 2;
                
                if(! m_A_sup.isOnBoundary(index3)){
                    vec[index3] = 0.5 * (vec[index1] + vec[index2]);
                }else{
                    continue;
                }
            }
            
            size_t width = m_A_sup.getWidth();
            for(size_t i = 0; i < m_A_sup.rows() / width; i++){
                for(size_t j = i * width; j < (i + 1) * width - 1; j += 2){
                    if(! m_A_sup.isOnBoundary(m_A_sup.mask(j+1)))
                        vec[m_A_sup.mask(j + 1)] = 0.5 * (vec[m_A_sup.mask(j)] + vec[m_A_sup.mask(j + 2)]);
                    else
                        continue;
                }
            }
        }

        friend std::vector<double>& operator*(std::vector<double> &x_k, InterpolationClass &B)
        {
            B.interpolate(x_k);
            return x_k;
        }
    

};

/*
template<class Vector, class Smoother>
class SawtoothMGIteration{
    private:
        std::vector<PoissonMatrix<double>> &A_level;
        Vector &b;

        std::vector<double> res;
        std::vector<double> err;

        std::vector<std::unique_ptr<SmootherClass<std::vector<double>>>> iterations;
        std::vector<std::unique_ptr<InterpolationClass>> interpolators;

        std::unique_ptr<Residual<Vector>> RES;
        std::unique_ptr<Residual<std::vector<double>>> COARSE_RES;

        std::unique_ptr<Solver<std::vector<double>>> COARSE_SOLVER;
        int nu = 2;

        //I need this for the gif
        unsigned int counter = 0;
        std::vector<double> temp;
        std::string name;

        //std::vector<int> &nu;

    public:
        SawtoothMGIteration(std::vector<PoissonMatrix<double>> &matrices, Vector &knownVec): A_level(matrices), b(knownVec) {
            res = std::vector<double>(b.size(),0.);
            err = std::vector<double>(b.size(),0.);

            
            for(size_t j = 0; j < A_level.size(); j++){
                iterations.push_back(std::make_unique<Smoother>(A_level[j],res));
            }
            
            
            for(size_t j = 0; j < A_level.size() - 1; j++){
                interpolators.push_back(std::make_unique<InterpolationClass>(A_level[j+1],A_level[j]));
            }
            

            RES = std::make_unique<Residual<Vector>>(A_level.at(0),b,res);
            COARSE_RES = std::make_unique<Residual<std::vector<double>>>(A_level.back(),res);

            COARSE_SOLVER = std::make_unique<AMG::Solver<std::vector<double>>>((*iterations.back()),(*COARSE_RES),2000,0.5,1);

            //for the gif
            temp = std::vector<double>(b.size());
        }


        
        std::vector<double> formatVector(std::vector<double> &in, PoissonMatrix<double> &A){
            std::vector<double> temp;
            for(size_t i = 0; i < A.rows(); i++){
                temp.push_back(in[A.mask(i)]);
            }
            return temp;
        }
        

        void apply_iteration_to_vec(std::vector<double> &sol) {
            //for the gif
            temp = sol;
            //********

            sol * (*RES);
            COARSE_RES->refresh_normalization_constant();

            //gif
            name = "./output/" + std::to_string(counter) + ".mtx";
            Utils::saveVectorOnFile(formatVector(temp,A_level.back()),name);
            counter++;
            //******

            err * (*COARSE_SOLVER) * (*COARSE_RES);
            std::cout<<"Achieved residual on coarse grid: "<<COARSE_RES->Norm()<<std::endl;

            //gif
            name = "./output/" + std::to_string(counter) + ".mtx";
            for(size_t j = 0; j < A_level.back().rows(); j++){
                auto id = A_level.back().mask(j);
                temp[id] += err[id];
            }
            Utils::saveVectorOnFile(formatVector(temp,A_level.back()),name);
            counter++;
            temp = sol;
            //******

            for(size_t j = A_level.size() - 1; j > 0; --j){
                err * (*interpolators[j-1]);
                //gif
                
                name = "./output/" + std::to_string(counter) + ".mtx";
                
                for(size_t k = 0; k < A_level.at(j-1).rows(); k++){
                    auto id = A_level.at(j-1).mask(k);
                    temp[id] += err[id];
                }
                

                Utils::saveVectorOnFile(formatVector(temp,A_level.at(j-1)),name);
                counter++;
                temp = sol;
                
                //******

                for(int i = 0; i < nu; i++){
                    err * (*iterations[j-1]);
                }

                //gif
                name = "./output/" + std::to_string(counter) + ".mtx";
                
                for(size_t k = 0; k < A_level.at(j-1).rows(); k++){
                    auto id = A_level.at(j-1).mask(k);
                    temp[id] += err[id];
                }
                

                Utils::saveVectorOnFile(formatVector(temp,A_level.at(j-1)),name);
                counter++;
                temp = sol;
                //******
            }
            
            for(size_t j = 0; j < sol.size(); j++){
                sol[j] += err[j];
                err[j] = 0;
            }
            
            //gif
            name = "./output/" + std::to_string(counter) + ".mtx";
            Utils::saveVectorOnFile(formatVector(sol,A_level.at(0)),name);
            counter++;
            //******

            
        }

        friend std::vector<double>& operator*(std::vector<double> &x_k, SawtoothMGIteration &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }



        ~SawtoothMGIteration(){
        }
};
*/


template<class Vector, class Smoother>
class SawtoothMGIteration{
    private:
        std::vector<PoissonMatrix<double>> &A_level;
        Vector &b;

        std::vector<double> res;
        std::vector<double> err;

        std::vector<std::unique_ptr<SmootherClass<std::vector<double>>>> iterations;
        std::vector<std::unique_ptr<InterpolationClass>> interpolators;

        std::unique_ptr<Residual<Vector>> RES;
        std::unique_ptr<Residual<std::vector<double>>> COARSE_RES;

        std::unique_ptr<Solver<std::vector<double>>> COARSE_SOLVER;
        int nu = 5;
        //std::vector<int> &nu;

    public:
        SawtoothMGIteration(std::vector<PoissonMatrix<double>> &matrices, Vector &knownVec): A_level(matrices), b(knownVec) {
            res = std::vector<double>(b.size(),0.);
            err = std::vector<double>(b.size(),0.);
          
            for(size_t j = 0; j < A_level.size(); j++){
                iterations.push_back(std::make_unique<Smoother>(A_level[j],res));
            }
                  
            for(size_t j = 0; j < A_level.size() - 1; j++){
                interpolators.push_back(std::make_unique<InterpolationClass>(A_level[j+1],A_level[j]));
            }    

            RES = std::make_unique<Residual<Vector>>(A_level.at(0),b,res);
            COARSE_RES = std::make_unique<Residual<std::vector<double>>>(A_level.back(),res);

            COARSE_SOLVER = std::make_unique<AMG::Solver<std::vector<double>>>((*iterations.back()),(*COARSE_RES),2000,1.e-1,1);
        }        

        void apply_iteration_to_vec(std::vector<double> &sol) {
            sol * (*RES);
            COARSE_RES->refresh_normalization_constant();

            err * (*COARSE_SOLVER) * (*COARSE_RES);
            std::cout<<"Achieved residual on coarse grid: "<<COARSE_RES->Norm()<<std::endl;


            for(size_t j = A_level.size() - 1; j > 0; --j){
                err * (*interpolators[j-1]);
                for(int i = 0; i < nu; i++){
                    err * (*iterations[j-1]);
                }
            }
            
            for(size_t j = 0; j < sol.size(); j++){
                sol[j] += err[j];
                err[j] = 0;
            }

            
        }

        friend std::vector<double>& operator*(std::vector<double> &x_k, SawtoothMGIteration &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }



        ~SawtoothMGIteration(){
        }
};


}


#endif

