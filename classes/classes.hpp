#ifndef CLASSES_H
#define CLASSES_H


/*
The main idea is to create an abstract class iteration, which takes as constructor a pointer to the matrix A, to the vector b, some data to
identify the current level grid nodes, to apply the iteration to x^k we could define an operator i.e.:
regularCoarsening = true;
GaussSeidelIteration GS(A,b,regularCoarsening,"h");
x-(GS)->-(GS)->
in this case for example we applied two gauss seidel iterations to an initial guess
*/

template<typename T_m, typename T_v>
class Iteration{
    public:
        virtual Iteration(){

        }
    protected:
        T_m* m_A;
        T_v* m_b;
        size_t m_size;
};

#endif