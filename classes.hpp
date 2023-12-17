#ifndef CLASSES_H
#define CLASSES_H

#include <tuple>
#include <vector>
#include <functional>

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
        virtual const std::vector<size_t> inRowConnections(const size_t l) const = 0;

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
            for(size_t i = 0; i < level; i++){
                width = (width + 1) / 2;
                step *= 2;
            }
        }

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

        const std::vector<size_t> inRowConnections(const size_t l) const override{
            auto equivalent_l = mask(l);
            std::vector<size_t> temp;
			if(isOnBoundary(equivalent_l)){
				temp.push_back(l);
			}else{
				temp.push_back(l - width);
				temp.push_back(l - 1);
				temp.push_back(l);
				temp.push_back(l + 1);
				temp.push_back(l + width);
			}
			return temp;
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
		
		const std::vector<size_t> nonZerosInRow(const size_t row){
			return m_domain.inRowConnections(row);
		}

        const size_t nonZeros(){
            return m_size + m_domain.numConnections();
        }

        const size_t mask(const size_t l){
            return m_domain.mask(l);
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

}


#endif